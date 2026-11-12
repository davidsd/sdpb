// Resident mode: sdpb --controlDir <dir>. See Readme.md.

#include "Resident.hxx"
#include "Control_Channel.hxx"
#include "sdpb/create_block_info.hxx"
#include "sdpb/save_solution.hxx"

#include "pmp/max_normalization_index.hxx"
#include "pmp2sdp/Block_File_Format.hxx"
#include "pmp2sdp/Dual_Constraint_Group.hxx"
#include "pmp2sdp/Output_SDP/Output_SDP.hxx"
#include "pmp2sdp/Reuse_SDP.hxx"
#include "pmp_read/PMP_File_Parse_Result.hxx"
#include "pmp_read/pmp_read.hxx"
#include "sdp_solve/SDP/read_block_data/Block_Data_Parse_Result.hxx"
#include "sdp_solve/SDP/read_block_data/SDP_Block_Data.hxx"
#include "sdp_solve/SDP/set_bases_blocks.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <csignal>
#include <cstdlib>
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <unistd.h>

#include <chrono>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <sstream>

namespace fs = std::filesystem;

// Defined elsewhere in pmp2sdp_lib / sdp_solve
void write_block_data(std::ostream &os, const Dual_Constraint_Group &group,
                      Block_File_Format format);
void write_objectives_json(std::ostream &output_stream,
                           const El::BigFloat &objective_const,
                           const std::vector<El::BigFloat> &dual_objective_b);
void write_normalization_json(std::ostream &output_stream,
                              const std::vector<El::BigFloat> &normalization);
Block_Data_Parse_Result parse_block_data(std::istream &block_stream,
                                         Block_File_Format format);
void read_objectives_from_stream(
  const El::Grid &grid, std::istream &objectives_stream,
  El::BigFloat &objective_const,
  El::DistMatrix<El::BigFloat> &dual_objective_b);
void read_objectives(const fs::path &sdp_path, const El::Grid &grid,
                     El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b,
                     Timers &timers);
std::optional<std::vector<El::BigFloat>>
read_normalization(const fs::path &sdp_path, Timers &timers);

namespace
{
  using Clock = std::chrono::high_resolution_clock;
  using Json_Writer = rapidjson::Writer<rapidjson::StringBuffer>;

  double seconds_since(const Clock::time_point &start)
  {
    return std::chrono::duration<double>(Clock::now() - start).count();
  }

  std::string read_file(const fs::path &path)
  {
    std::ifstream in(path, std::ios::binary);
    ASSERT(in.good(), "Cannot open ", path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  void broadcast_string(std::string &s, const int root = 0,
                        const El::mpi::Comm &comm = El::mpi::COMM_WORLD)
  {
    size_t size = s.size();
    El::mpi::Broadcast(size, root, comm);
    s.resize(size);
    if(size > 0)
      MPI_Bcast(s.data(), size, MPI_CHAR, root, comm.comm);
  }

  // An error that all ranks see the same way (validation before any change
  // to the state): the command fails and SDPB carries on. Any other exception
  // aborts MPI, as in sdpb, since some ranks may be inside a collective.
  struct Command_Error : std::runtime_error
  {
    using std::runtime_error::runtime_error;
  };

  // Collective: if any rank has a non-empty error, all ranks throw the error
  // of the lowest such rank as a Command_Error.
  void agree_or_throw(const std::string &error)
  {
    int first = error.empty() ? El::mpi::Size() : El::mpi::Rank();
    first = El::mpi::AllReduce(first, El::mpi::MIN, El::mpi::COMM_WORLD);
    if(first == El::mpi::Size())
      return;
    std::string message = error;
    broadcast_string(message, first);
    throw Command_Error(message);
  }

  // Run f (which must not call collectives) on every rank, then agree on
  // whether it failed. Collective.
  template <class F> void consistently(F &&f)
  {
    std::string error;
    try
      {
        f();
      }
    catch(std::exception &e)
      {
        error = e.what();
        if(error.empty())
          error = "unknown error";
      }
    agree_or_throw(error);
  }

  std::string to_string(const El::BigFloat &x)
  {
    std::ostringstream os;
    set_stream_precision(os);
    os << x;
    return os.str();
  }

  // A decimal number, or a rational p/q
  El::BigFloat parse_number(const std::string &s)
  {
    const auto slash = s.find('/');
    if(slash == std::string::npos)
      return El::BigFloat(s);
    return El::BigFloat(s.substr(0, slash))
           / El::BigFloat(s.substr(slash + 1));
  }

  std::string number_string(const rapidjson::Value &v)
  {
    if(v.IsString())
      return v.GetString();
    ASSERT(v.IsNumber(), "Expected a number or a string");
    // Exact for integers; give other values as strings
    ASSERT(v.IsInt64(), "Give non-integer numbers as strings");
    return std::to_string(v.GetInt64());
  }

  std::string get_string(const rapidjson::Value &v, const char *key,
                         const std::string &default_value)
  {
    if(!v.HasMember(key))
      return default_value;
    ASSERT(v[key].IsString(), "'", key, "' must be a string");
    return v[key].GetString();
  }

  bool get_bool(const rapidjson::Value &v, const char *key,
                const bool default_value)
  {
    if(!v.HasMember(key))
      return default_value;
    ASSERT(v[key].IsBool(), "'", key, "' must be true or false");
    return v[key].GetBool();
  }

  std::vector<std::string>
  get_string_list(const rapidjson::Value &v, const char *key)
  {
    std::vector<std::string> result;
    if(!v.HasMember(key))
      return result;
    if(v[key].IsString())
      {
        result.emplace_back(v[key].GetString());
        return result;
      }
    ASSERT(v[key].IsArray(), "'", key, "' must be a string or an array");
    for(const auto &x : v[key].GetArray())
      result.push_back(number_string(x));
    return result;
  }

  void write_bigfloat(Json_Writer &w, const El::BigFloat &x)
  {
    const auto s = to_string(x);
    w.String(s.c_str(), s.size());
  }

  // Value of option name (or short_name) in a space-separated command line
  std::optional<std::string>
  find_option(const std::vector<std::string> &tokens, const std::string &name,
              const std::string &short_name)
  {
    for(size_t i = 0; i < tokens.size(); ++i)
      {
        const auto &t = tokens.at(i);
        if(t == name || t == short_name)
          {
            ASSERT(i + 1 < tokens.size(), "Missing value for ", t);
            return tokens.at(i + 1);
          }
        if(t.rfind(name + "=", 0) == 0)
          return t.substr(name.size() + 1);
      }
    return std::nullopt;
  }

  // What the resident process needs to know about the files of an SDP
  // directory written by pmp2sdp. Read on rank 0, the same on all ranks.
  struct Sdp_Files
  {
    fs::path dir;
    size_t num_blocks = 0;
    int64_t max_num_poles = -1;
    Block_File_Format format = Block_File_Format::bin;
    // Source PMP file of each block (pmp_info.json), empty if unknown
    std::vector<fs::path> block_paths;
    // Text of objectives.json and normalization.json ("" if absent)
    std::string objectives_text, normalization_text;
    bool is_directory = true;

    explicit Sdp_Files(const fs::path &sdp_dir) : dir(sdp_dir)
    {
      std::string control, pmp_info;
      int dir_flag = 1;
      std::string error;
      if(El::mpi::Rank() == 0)
        {
          try
            {
              dir_flag = fs::is_directory(sdp_dir) ? 1 : 0;
              if(dir_flag)
                {
                  control = read_file(sdp_dir / "control.json");
                  if(fs::exists(sdp_dir / "pmp_info.json"))
                    pmp_info = read_file(sdp_dir / "pmp_info.json");
                  objectives_text = read_file(sdp_dir / "objectives.json");
                  if(fs::exists(sdp_dir / "normalization.json"))
                    normalization_text
                      = read_file(sdp_dir / "normalization.json");
                }
            }
          catch(std::exception &e)
            {
              error = e.what();
            }
        }
      agree_or_throw(error);
      El::mpi::Broadcast(dir_flag, 0, El::mpi::COMM_WORLD);
      is_directory = dir_flag != 0;
      if(!is_directory)
        return;
      broadcast_string(control);
      broadcast_string(pmp_info);
      broadcast_string(objectives_text);
      broadcast_string(normalization_text);

      // The same text on all ranks: errors are the same everywhere
      consistently([&] { parse(control, pmp_info, sdp_dir); });
    }

  private:
    void parse(const std::string &control, const std::string &pmp_info,
               const fs::path &sdp_dir)
    {
      rapidjson::Document control_doc;
      control_doc.Parse(control.c_str());
      ASSERT(!control_doc.HasParseError() && control_doc.IsObject()
               && control_doc.HasMember("num_blocks")
               && control_doc["num_blocks"].IsUint64(),
             "Cannot parse ", sdp_dir / "control.json");
      num_blocks = control_doc["num_blocks"].GetUint64();
      if(control_doc.HasMember("command") && control_doc["command"].IsString())
        {
          std::vector<std::string> tokens;
          std::istringstream iss(control_doc["command"].GetString());
          std::string token;
          while(iss >> token)
            tokens.push_back(token);
          const auto poles
            = find_pmp2sdp_option(tokens, "--maxNumPoles", "-n");
          if(poles.has_value())
            max_num_poles = std::stoll(poles.value());
          const auto output_format
            = find_pmp2sdp_option(tokens, "--outputFormat", "-f");
          if(output_format.has_value() && output_format.value() == "json")
            format = Block_File_Format::json;
        }
      if(!pmp_info.empty())
        {
          rapidjson::Document doc;
          doc.Parse(pmp_info.c_str());
          ASSERT(!doc.HasParseError() && doc.IsArray(), "Cannot parse ",
                 sdp_dir / "pmp_info.json");
          block_paths.resize(doc.Size());
          for(const auto &entry : doc.GetArray())
            {
              ASSERT(entry.IsObject() && entry.HasMember("index")
                       && entry["index"].IsUint64()
                       && entry["index"].GetUint64() < block_paths.size()
                       && entry.HasMember("path") && entry["path"].IsString(),
                     "Cannot parse ", sdp_dir / "pmp_info.json");
              block_paths.at(entry["index"].GetUint64())
                = fs::path(entry["path"].GetString());
            }
          ASSERT_EQUAL(block_paths.size(), num_blocks,
                       sdp_dir / "pmp_info.json");
        }
    }
  };

  struct Iterate
  {
    Block_Vector x, y;
    Paired_Block_Diagonal_Matrix X, Y;
  };

  // A block whose data is c = sum_i w_i c_i, B = sum_i w_i B_i
  struct Family
  {
    size_t block_index = 0;
    int owner_root = 0;
    size_t num_components = 0;
    // A symmetric form_dim x form_dim matrix of components, given by its
    // upper triangle (row-major), if form_dim > 0
    size_t form_dim = 0;
    std::vector<std::string> paths;
    // Only on owner_root:
    std::vector<El::Matrix<El::BigFloat>> c, B;
    std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases;
    // Sample scaling of the first sample point of each component
    std::vector<El::BigFloat> first_scaling;
    // normalization.json text of the normalization the components were
    // converted with (set_coefficients refuses another)
    std::string normalization_text;
  };

  class Resident_State
  {
  public:
    Resident_State(const Environment &env, const SDPB_Parameters &parameters,
                   const std::vector<std::string> &startup_args,
                   mp_bitcnt_t startup_precision,
                   const Block_Info &block_info);

    std::string ready_json(double load_seconds) const;
    // Execute one command and write its result (an object). Returns false if
    // the resident loop should stop.
    bool execute(const rapidjson::Value &command, Json_Writer &w);
    void checkpoint(const fs::path &dir,
                    const boost::property_tree::ptree &properties);
    const SDPB_Parameters &startup_parameters() const { return parameters; }
    // The iterate in memory is newer than the last checkpoint written
    bool dirty = false;
    // A solve stopped on SIGTERM
    bool sigterm = false;

  private:
    const Environment &env;
    SDPB_Parameters parameters;
    std::vector<std::string> startup_args;
    mp_bitcnt_t startup_precision;
    const Block_Info &block_info;
    El::Grid grid;
    Timers load_timers;
    SDP sdp;
    SDP_Solver solver;
    bool loaded_checkpoint = false;
    std::unique_ptr<BigInt_Shared_Memory_Syrk_Context> syrk_context;
    // World rank of the root of the group that owns each block
    std::vector<int> owner_root;
    Sdp_Files files;
    // Normalization, on all ranks (sdp.normalization is on rank 0 only)
    std::optional<std::vector<El::BigFloat>> normalization;
    fs::path block_timings_path;
    std::map<std::string, Iterate> iterates;
    std::map<std::string, Family> families;
    size_t num_solves = 0;

    bool is_group_root() const { return grid.Rank() == 0; }
    size_t num_blocks() const { return block_info.dimensions.size(); }
    std::optional<size_t> local_index(size_t block_index) const;

    void cmd_solve(const rapidjson::Value &command, Json_Writer &w,
                   bool &keep_running);
    void cmd_update_sdp(const rapidjson::Value &command, Json_Writer &w);
    void cmd_update_pmp(const rapidjson::Value &command, Json_Writer &w);
    void cmd_load_family(const rapidjson::Value &command, Json_Writer &w);
    void cmd_set_coefficients(const rapidjson::Value &command,
                              Json_Writer &w);
    void cmd_negate_normalization(Json_Writer &w);
    void cmd_save_iterate(const rapidjson::Value &command, Json_Writer &w);
    void cmd_restore_iterate(const rapidjson::Value &command, Json_Writer &w);
    void cmd_status(Json_Writer &w);

    std::map<size_t, std::string>
    check_new_blocks(const std::vector<Dual_Constraint_Group> &groups,
                     const std::vector<size_t> &changed,
                     std::vector<int> &holder) const;
    void replace_blocks(std::map<size_t, std::string> bytes,
                        const std::vector<size_t> &changed,
                        const std::vector<int> &holder);
    void check_checkpoint_dir(const fs::path &dir) const;
    void checkpoint_elsewhere(const fs::path &dir);
    void set_objectives_text(const std::string &text);
    El::Matrix<El::BigFloat> y_on_group_root() const;
    void write_derived(const std::vector<std::string> &names,
                       Json_Writer &w) const;
  };

  Resident_State::Resident_State(const Environment &env,
                                 const SDPB_Parameters &parameters,
                                 const std::vector<std::string> &startup_args,
                                 const mp_bitcnt_t startup_precision,
                                 const Block_Info &block_info)
      : env(env),
        parameters(parameters),
        startup_args(startup_args),
        startup_precision(startup_precision),
        block_info(block_info),
        grid(block_info.mpi_comm.value),
        load_timers(env, parameters.verbosity),
        sdp(parameters.sdp_path, block_info, grid, load_timers,
            parameters.solver.sparse_B),
        solver(parameters.solver, parameters.verbosity,
               parameters.require_initial_checkpoint, block_info, grid,
               sdp.dual_objective_b.Height()),
        files(parameters.sdp_path)
  {
    loaded_checkpoint = solver.current_generation != 0
                        || solver.backup_generation.has_value();
    owner_root.assign(num_blocks(), -1);
    if(is_group_root())
      for(const auto index : block_info.block_indices)
        owner_root.at(index) = El::mpi::Rank();
    El::mpi::AllReduce(owner_root.data(), owner_root.size(), El::mpi::MAX,
                       El::mpi::COMM_WORLD);

    // Normalization on all ranks
    {
      std::vector<El::BigFloat> n;
      int has = 0;
      if(El::mpi::Rank() == 0 && sdp.normalization.has_value())
        {
          n = sdp.normalization.value();
          has = 1;
        }
      El::mpi::Broadcast(has, 0, El::mpi::COMM_WORLD);
      if(has)
        {
          size_t size = n.size();
          El::mpi::Broadcast(size, 0, El::mpi::COMM_WORLD);
          n.resize(size);
          El::mpi::Broadcast(n.data(), size, 0, El::mpi::COMM_WORLD);
          normalization = n;
        }
    }
    block_timings_path = get_block_timings_path(parameters);
    if(block_timings_path.empty()
       && fs::exists(parameters.solver.checkpoint_out / "block_timings"))
      block_timings_path = parameters.solver.checkpoint_out / "block_timings";
  }

  std::optional<size_t> Resident_State::local_index(size_t block_index) const
  {
    const auto &indices = block_info.block_indices;
    const auto it = std::find(indices.begin(), indices.end(), block_index);
    if(it == indices.end())
      return std::nullopt;
    return std::distance(indices.begin(), it);
  }

  std::string Resident_State::ready_json(const double load_seconds) const
  {
    rapidjson::StringBuffer sb;
    Json_Writer w(sb);
    w.StartObject();
    w.Key("status");
    w.String("ready");
    char host[256] = {0};
    gethostname(host, sizeof(host) - 1);
    w.Key("host");
    w.String(host);
    w.Key("pid");
    w.Int64(getpid());
    w.Key("sdpDir");
    w.String(parameters.sdp_path.string().c_str());
    w.Key("num_blocks");
    w.Uint64(num_blocks());
    w.Key("dual_dimension");
    w.Int64(sdp.dual_objective_b.Height());
    w.Key("precision");
    w.Uint64(parameters.solver.precision);
    w.Key("mpi_processes");
    w.Int(El::mpi::Size());
    w.Key("nodes");
    w.Int(env.num_nodes());
    w.Key("loaded_checkpoint");
    w.Bool(loaded_checkpoint);
    w.Key("load_seconds");
    w.Double(load_seconds);
    // Protocol features a driver may check for
    w.Key("features");
    w.StartArray();
    w.String("solve_inherit");
    w.EndArray();
    w.EndObject();
    return sb.GetString();
  }

  void
  Resident_State::checkpoint(const fs::path &dir,
                             const boost::property_tree::ptree &properties)
  {
    solver.save_checkpoint(dir, parameters.verbosity, properties);
    if(dir == parameters.solver.checkpoint_out)
      dirty = false;
    // Keep the block distribution with the checkpoint, as the file route
    // does, so that a restart from it gives the same layout
    if(El::mpi::Rank() == 0 && !block_timings_path.empty()
       && fs::exists(block_timings_path)
       && !fs::exists(dir / "block_timings"))
      fs::copy_file(block_timings_path, dir / "block_timings");
  }

  // A checkpoint in a directory other than --checkpointDir: keep the
  // generation counters of the main one
  void Resident_State::checkpoint_elsewhere(const fs::path &dir)
  {
    const auto current = solver.current_generation;
    const auto backup = solver.backup_generation;
    solver.current_generation = 0;
    solver.backup_generation.reset();
    checkpoint(dir, to_property_tree(parameters));
    solver.current_generation = current;
    solver.backup_generation = backup;
  }

  // Collective: Command_Error unless dir holds a binary checkpoint with a
  // file for every rank (so that loading it does not fail on some ranks)
  void Resident_State::check_checkpoint_dir(const fs::path &dir) const
  {
    int64_t generation = -1;
    std::string error;
    if(El::mpi::Rank() == 0)
      {
        try
          {
            const auto metadata = dir / "checkpoint.json";
            if(fs::exists(metadata))
              {
                rapidjson::Document doc;
                doc.Parse(read_file(metadata).c_str());
                ASSERT(!doc.HasParseError() && doc.IsObject()
                         && doc.HasMember("current")
                         && doc["current"].IsInt64(),
                       "Cannot parse ", metadata);
                generation = doc["current"].GetInt64();
              }
          }
        catch(std::exception &e)
          {
            error = e.what();
          }
      }
    agree_or_throw(error);
    El::mpi::Broadcast(generation, 0, El::mpi::COMM_WORLD);
    const auto file
      = generation >= 0
          ? dir
              / ("checkpoint_" + std::to_string(generation) + "_"
                 + std::to_string(El::mpi::Rank()))
          : dir / ("checkpoint." + std::to_string(El::mpi::Rank()));
    agree_or_throw(fs::exists(file)
                     ? std::string()
                     : "No checkpoint file " + file.string()
                         + " (a checkpoint for this many ranks?)");
  }

  El::Matrix<El::BigFloat> Resident_State::y_on_group_root() const
  {
    // y is duplicated among the blocks: gather the first local block's copy
    // to the root of each group (collective in each group)
    const auto &y_dist = solver.y.blocks.at(0);
    El::DistMatrix<El::BigFloat, El::CIRC, El::CIRC> y_circ(y_dist);
    if(is_group_root())
      return y_circ.Matrix();
    return {};
  }

  bool Resident_State::execute(const rapidjson::Value &command, Json_Writer &w)
  {
    // The same text on all ranks: these errors are the same everywhere
    std::string name;
    consistently([&] {
      ASSERT(command.IsObject(), "A command must be a JSON object");
      name = get_string(command, "cmd", "");
      ASSERT(name == "solve" || name == "update_sdp" || name == "update_pmp"
               || name == "load_family" || name == "set_coefficients"
               || name == "negate_normalization" || name == "save_iterate"
               || name == "restore_iterate" || name == "drop_iterate"
               || name == "checkpoint" || name == "status" || name == "quit",
             "Unknown command: '", name, "'");
    });
    bool keep_running = true;
    const auto start = Clock::now();
    w.StartObject();
    w.Key("cmd");
    w.String(name.c_str());
    if(name == "solve")
      cmd_solve(command, w, keep_running);
    else if(name == "update_sdp")
      cmd_update_sdp(command, w);
    else if(name == "update_pmp")
      cmd_update_pmp(command, w);
    else if(name == "load_family")
      cmd_load_family(command, w);
    else if(name == "set_coefficients")
      cmd_set_coefficients(command, w);
    else if(name == "negate_normalization")
      cmd_negate_normalization(w);
    else if(name == "save_iterate")
      cmd_save_iterate(command, w);
    else if(name == "restore_iterate")
      cmd_restore_iterate(command, w);
    else if(name == "drop_iterate")
      {
        std::string iterate_name;
        consistently([&] { iterate_name = get_string(command, "name", ""); });
        iterates.erase(iterate_name);
      }
    else if(name == "checkpoint")
      {
        fs::path dir;
        consistently([&] {
          dir = get_string(command, "dir",
                           parameters.solver.checkpoint_out.string());
        });
        if(dir == parameters.solver.checkpoint_out)
          checkpoint(dir, to_property_tree(parameters));
        else
          checkpoint_elsewhere(dir);
        w.Key("dir");
        w.String(dir.string().c_str());
      }
    else if(name == "status")
      cmd_status(w);
    else if(name == "quit")
      {
        bool with_checkpoint = false;
        consistently(
          [&] { with_checkpoint = get_bool(command, "checkpoint", false); });
        if(with_checkpoint)
          checkpoint(parameters.solver.checkpoint_out,
                     to_property_tree(parameters));
        keep_running = false;
      }
    El::mpi::Barrier(El::mpi::COMM_WORLD);
    w.Key("seconds");
    w.Double(seconds_since(start));
    w.EndObject();
    return keep_running;
  }

  // solve: run the solver from the current (or a named) iterate.
  // "args": an sdpb command line for this solve (outDir, checkpointDir,
  // thresholds, maxIterations, ...); options it lacks come from the resident
  // process's command line, unless "inherit" is false: then the solve's
  // options are its own args only, exactly as for an sdpb run with that
  // command line (the options fixed at start-up are still checked).
  void Resident_State::cmd_solve(const rapidjson::Value &command,
                                 Json_Writer &w, bool &keep_running)
  {
    const auto start = Clock::now();
    // Validate everything before the solver runs, the same way on all ranks
    // (errors on rank 0 only, e.g. creating --outDir, are shared)
    std::vector<std::string> args, derived;
    std::optional<SDPB_Parameters> parsed;
    std::string start_kind = "current", start_iterate;
    bool with_checkpoint = false, with_solution = true, inherit = true;
    consistently([&] {
      args = get_string_list(command, "args");
      derived = get_string_list(command, "derived");
      with_checkpoint = get_bool(command, "checkpoint", false);
      with_solution = get_bool(command, "write_solution", true);
      inherit = get_bool(command, "inherit", true);
      // The SDP is the one in memory: a different --sdpDir is a mistake
      // (decided on rank 0: file system checks may differ between ranks)
      const auto sdp_dir = find_option(args, "--sdpDir", "-s");
      if(sdp_dir.has_value() && El::mpi::Rank() == 0)
        {
          std::error_code ec;
          const bool same
            = fs::path(sdp_dir.value()).lexically_normal()
                == files.dir.lexically_normal()
              || fs::equivalent(sdp_dir.value(), files.dir, ec);
          ASSERT(same, "solve: --sdpDir ", sdp_dir.value(),
                 " is not the SDP in memory (", files.dir,
                 "); use update_sdp or update_pmp to change it");
        }
      // Parse at the precision sdpb parses its command line at
      const mp_bitcnt_t working_precision = mpf_get_default_prec();
      mpf_set_default_prec(startup_precision);
      try
        {
          if(inherit)
            parsed.emplace(
              std::vector<std::vector<std::string>>{args, startup_args},
              true, env);
          else
            parsed.emplace(std::vector<std::vector<std::string>>{args}, true,
                           env);
        }
      catch(...)
        {
          mpf_set_default_prec(working_precision);
          throw;
        }
      mpf_set_default_prec(working_precision);
      const SDPB_Parameters &q = parsed.value();
      ASSERT_EQUAL(q.solver.precision, parameters.solver.precision,
                   "solve: --precision is fixed in resident mode");
      ASSERT_EQUAL(q.solver.sparse_B, parameters.solver.sparse_B,
                   "solve: --sparseB is fixed in resident mode");
      // The BigInt syrk context, sized by these limits, is kept between
      // solves
      ASSERT_EQUAL(q.solver.max_shared_memory.bytes,
                   parameters.solver.max_shared_memory.bytes,
                   "solve: --maxSharedMemory is fixed in resident mode");
      ASSERT_EQUAL(q.solver.max_memory.bytes, parameters.solver.max_memory.bytes,
                   "solve: --maxMemory is fixed in resident mode");
      ASSERT_EQUAL(q.proc_granularity, parameters.proc_granularity,
                   "solve: --procGranularity is fixed in resident mode");
      // The kept syrk context is sized with the factor's memory
      ASSERT_EQUAL(q.solver.multimodular_trsm,
                   parameters.solver.multimodular_trsm,
                   "solve: --multimodularTrsm is fixed in resident mode");
      // Starting point
      if(command.HasMember("start"))
        {
          const auto &v = command["start"];
          if(v.IsString())
            start_kind = v.GetString();
          else
            {
              ASSERT(v.IsObject(),
                     "solve: 'start' must be a string or {\"iterate\": name}");
              start_kind = "iterate";
              start_iterate = get_string(v, "iterate", "");
              ASSERT(iterates.count(start_iterate) != 0,
                     "solve: no iterate named '", start_iterate, "'");
            }
          ASSERT(start_kind == "current" || start_kind == "initial"
                   || start_kind == "checkpoint" || start_kind == "iterate",
                 "solve: unknown start '", start_kind, "'");
        }
      for(const auto &name : derived)
        ASSERT(name == "y" || name == "z" || name == "family_values",
               "solve: unknown derived quantity '", name, "'");
      const bool needs_z
        = std::count(derived.begin(), derived.end(), "z") != 0
          || (with_solution && q.write_solution.vector_z);
      ASSERT(!needs_z || normalization.has_value(),
             "solve: z needs a normalization");
    });
    const SDPB_Parameters &p = parsed.value();
    if(start_kind == "checkpoint")
      check_checkpoint_dir(p.solver.checkpoint_in);

    // Starting point (validated above)
    if(start_kind == "initial")
      {
        solver.X.set_zero();
        solver.Y.set_zero();
        for(auto &block : solver.x.blocks)
          El::Zero(block);
        for(auto &block : solver.y.blocks)
          El::Zero(block);
        solver.X.add_diagonal(p.solver.initial_matrix_scale_primal);
        solver.Y.add_diagonal(p.solver.initial_matrix_scale_dual);
        dirty = true;
      }
    else if(start_kind == "checkpoint")
      {
        solver.load_checkpoint(p.solver.checkpoint_in, block_info, p.verbosity,
                               true);
        dirty = true;
      }
    else if(start_kind == "iterate")
      {
        const auto &it = iterates.at(start_iterate);
        solver.x = it.x;
        solver.X = it.X;
        solver.y = it.y;
        solver.Y = it.Y;
        dirty = true;
      }

    if(p.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
      El::Output(boost::posix_time::second_clock::local_time(),
                 " Resident solve ", num_solves + 1, ", outDir ",
                 p.out_directory);

    Timers timers(env, p.verbosity);
    const boost::property_tree::ptree parameters_tree(to_property_tree(p));
    const auto iterations_json_path = p.out_directory / "iterations.json";
    El::Matrix<int32_t> block_timings_ms;
    const auto run_start = Clock::now();
    const SDP_Solver_Terminate_Reason reason(solver.run(
      env, p.solver, p.verbosity, parameters_tree, block_info, sdp, grid,
      run_start, iterations_json_path, timers, block_timings_ms,
      &syrk_context));
    const double run_seconds = seconds_since(run_start);
    ++num_solves;
    dirty = true;

    if(p.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
      {
        set_stream_precision(std::cout);
        std::cout << "-----" << reason << "-----\n"
                  << '\n'
                  << "primalObjective = " << solver.primal_objective << '\n'
                  << "dualObjective   = " << solver.dual_objective << '\n'
                  << "dualityGap      = " << solver.duality_gap << '\n'
                  << "primalError     = " << solver.primal_error() << '\n'
                  << "dualError       = " << solver.dual_error << '\n'
                  << '\n';
      }

    const bool sigterm_now
      = reason == SDP_Solver_Terminate_Reason::SIGTERM_Received;
    const auto checkpoint_start = Clock::now();
    if(sigterm_now || with_checkpoint)
      checkpoint(p.solver.checkpoint_out, parameters_tree);
    const double checkpoint_seconds = seconds_since(checkpoint_start);

    const auto solution_start = Clock::now();
    if(with_solution)
      {
        const auto runtime = static_cast<int64_t>(seconds_since(run_start));
        save_solution(solver, block_info, sdp, p, reason, runtime);
      }
    const double solution_seconds = seconds_since(solution_start);

    std::ostringstream reason_ss;
    reason_ss << reason;
    w.Key("terminateReason");
    w.String(reason_ss.str().c_str());
    w.Key("primalObjective");
    write_bigfloat(w, solver.primal_objective);
    w.Key("dualObjective");
    write_bigfloat(w, solver.dual_objective);
    w.Key("dualityGap");
    write_bigfloat(w, solver.duality_gap);
    w.Key("primalError");
    write_bigfloat(w, solver.primal_error());
    w.Key("dualError");
    write_bigfloat(w, solver.dual_error);
    w.Key("steps");
    w.Int64(solver.num_steps);
    w.Key("outDir");
    w.String(p.out_directory.string().c_str());
    w.Key("solve_index");
    w.Uint64(num_solves);
    w.Key("timing");
    w.StartObject();
    w.Key("setup");
    w.Double(std::chrono::duration<double>(run_start - start).count());
    w.Key("run");
    w.Double(run_seconds);
    w.Key("checkpoint");
    w.Double(checkpoint_seconds);
    w.Key("solution");
    w.Double(solution_seconds);
    w.EndObject();
    write_derived(derived, w);

    if(sigterm_now)
      {
        sigterm = true;
        keep_running = false;
      }
  }

  // Derived quantities after a solve:
  //   "y", "z": the dual vector y, and z (y with the normalization restored);
  //   "family_values": for each family and component i, c_i - B_i.y over the
  //     block's constraint rows (for a 1x1 constant block: the functional on
  //     the component, times the positive sample scaling; "pmp_values"
  //     divides by it).
  // More can be added here (plug-in hook).
  void Resident_State::write_derived(const std::vector<std::string> &names,
                                     Json_Writer &w) const
  {
    if(names.empty())
      return;
    const auto y = y_on_group_root();
    w.Key("derived");
    w.StartObject();
    for(const auto &name : names)
      {
        if(name == "y" || name == "z")
          {
            // rank 0 is the root of its group, so it has y
            std::vector<std::string> values;
            if(El::mpi::Rank() == 0)
              {
                if(name == "y")
                  {
                    for(int i = 0; i < y.Height(); ++i)
                      values.push_back(to_string(y(i, 0)));
                  }
                else
                  {
                    ASSERT(normalization.has_value(),
                           "z needs a normalization");
                    const auto &n = normalization.value();
                    const auto max_index = max_normalization_index(n);
                    std::vector<El::BigFloat> z(y.Height() + 1);
                    for(size_t i = 0; i < max_index; ++i)
                      z.at(i) = y(i, 0);
                    z.at(max_index) = 0;
                    for(int i = max_index; i < y.Height(); ++i)
                      z.at(i + 1) = y(i, 0);
                    // Same arithmetic as save_solution()
                    El::Matrix<El::BigFloat> z_view, n_view;
                    z_view.LockedAttach(z.size(), 1, z.data(), z.size());
                    n_view.LockedAttach(n.size(), 1, n.data(), n.size());
                    auto nz = El::Dot(n_view, z_view);
                    z.at(max_index) = (1 - nz) / n.at(max_index);
                    for(const auto &x : z)
                      values.push_back(to_string(x));
                  }
                w.Key(name.c_str());
                w.StartArray();
                for(const auto &s : values)
                  w.String(s.c_str(), s.size());
                w.EndArray();
              }
          }
        else if(name == "family_values")
          {
            if(El::mpi::Rank() == 0)
              {
                w.Key("family_values");
                w.StartObject();
              }
            for(const auto &[family_name, family] : families)
              {
                // On the owner root: one JSON text per family
                std::string text;
                if(El::mpi::Rank() == family.owner_root)
                  {
                    rapidjson::StringBuffer sb;
                    Json_Writer fw(sb);
                    fw.StartObject();
                    fw.Key("block");
                    fw.Uint64(family.block_index);
                    std::vector<std::vector<El::BigFloat>> values;
                    for(size_t i = 0; i < family.num_components; ++i)
                      {
                        const auto &c = family.c.at(i);
                        const auto &B = family.B.at(i);
                        std::vector<El::BigFloat> v(c.Height());
                        for(int row = 0; row < c.Height(); ++row)
                          {
                            El::BigFloat By(0);
                            for(int col = 0; col < B.Width(); ++col)
                              By += B(row, col) * y(col, 0);
                            v.at(row) = c(row, 0) - By;
                          }
                        values.push_back(v);
                      }
                    fw.Key("values");
                    fw.StartArray();
                    for(const auto &v : values)
                      {
                        fw.StartArray();
                        for(const auto &x : v)
                          write_bigfloat(fw, x);
                        fw.EndArray();
                      }
                    fw.EndArray();
                    if(!values.empty() && values.front().size() == 1)
                      {
                        fw.Key("pmp_values");
                        fw.StartArray();
                        for(size_t i = 0; i < values.size(); ++i)
                          write_bigfloat(fw, values.at(i).at(0)
                                               / family.first_scaling.at(i));
                        fw.EndArray();
                      }
                    fw.EndObject();
                    text = sb.GetString();
                  }
                // Send to rank 0
                if(family.owner_root != 0)
                  {
                    if(El::mpi::Rank() == family.owner_root)
                      {
                        uint64_t size = text.size();
                        MPI_Send(&size, 1, MPI_UINT64_T, 0, 7,
                                 El::mpi::COMM_WORLD.comm);
                        MPI_Send(text.data(), size, MPI_CHAR, 0, 8,
                                 El::mpi::COMM_WORLD.comm);
                      }
                    else if(El::mpi::Rank() == 0)
                      {
                        uint64_t size = 0;
                        MPI_Recv(&size, 1, MPI_UINT64_T, family.owner_root, 7,
                                 El::mpi::COMM_WORLD.comm, MPI_STATUS_IGNORE);
                        text.resize(size);
                        MPI_Recv(text.data(), size, MPI_CHAR,
                                 family.owner_root, 8,
                                 El::mpi::COMM_WORLD.comm, MPI_STATUS_IGNORE);
                      }
                  }
                if(El::mpi::Rank() == 0)
                  {
                    w.Key(family_name.c_str());
                    w.RawValue(text.c_str(), text.size(),
                               rapidjson::kObjectType);
                  }
              }
            if(El::mpi::Rank() == 0)
              w.EndObject();
          }
        else
          RUNTIME_ERROR("Unknown derived quantity: '", name, "'");
      }
    w.EndObject();
  }

  // Collective. groups: the converted blocks held by this rank (any of them
  // may be unchanged); each changed block must be held by one rank. Serialize the changed blocks held by this rank and check that
  // every changed block is held by some rank and keeps its size (dim, number
  // of points, and N columns of B); Command_Error otherwise. holder: the rank
  // holding each block (-1: none).
  std::map<size_t, std::string> Resident_State::check_new_blocks(
    const std::vector<Dual_Constraint_Group> &groups,
    const std::vector<size_t> &changed, std::vector<int> &holder) const
  {
    const int rank = El::mpi::Rank();
    holder.assign(num_blocks(), -1);
    const std::set<size_t> changed_set(changed.begin(), changed.end());
    const auto N = sdp.dual_objective_b.Height();
    std::string error;
    std::map<size_t, std::string> bytes;
    for(const auto &group : groups)
      {
        if(changed_set.count(group.block_index) == 0)
          continue;
        holder.at(group.block_index) = rank;
        if(group.dim != block_info.dimensions.at(group.block_index)
           || group.num_points != block_info.num_points.at(group.block_index)
           || group.constraint_matrix.Width() != N)
          error = El::BuildString(
            "Block ", group.block_index,
            " would change its size (dim, points or N); the block structure "
            "must stay the same");
        std::ostringstream os(std::ios::binary);
        write_block_data(os, group, files.format);
        bytes[group.block_index] = os.str();
      }
    agree_or_throw(error);
    El::mpi::AllReduce(holder.data(), holder.size(), El::mpi::MAX,
                       El::mpi::COMM_WORLD);
    for(const auto index : changed)
      if(holder.at(index) < 0)
        throw Command_Error(El::BuildString("No data for block ", index));
    return bytes;
  }

  // Replace the changed blocks by the serialized blocks from
  // check_new_blocks(). The data go through the same serialization as a
  // block_data file (the SDP's --outputFormat), so the blocks in memory equal
  // those SDPB would read from pmp2sdp's output.
  void Resident_State::replace_blocks(std::map<size_t, std::string> bytes,
                                      const std::vector<size_t> &changed,
                                      const std::vector<int> &holder)
  {
    const int rank = El::mpi::Rank();
    const std::set<size_t> changed_set(changed.begin(), changed.end());
    // Send each block to the root of the group that owns it. Both sides go
    // through the blocks in increasing order, and MPI keeps the order of the
    // messages between two ranks, so two tags are enough.
    constexpr int size_tag = 1, data_tag = 2;
    std::vector<MPI_Request> requests;
    std::vector<uint64_t> sizes;
    sizes.reserve(bytes.size());
    for(const auto &[index, data] : bytes)
      {
        const int dest = owner_root.at(index);
        if(dest == rank)
          continue;
        sizes.push_back(data.size());
        requests.emplace_back();
        MPI_Isend(&sizes.back(), 1, MPI_UINT64_T, dest, size_tag,
                  El::mpi::COMM_WORLD.comm, &requests.back());
        requests.emplace_back();
        MPI_Isend(data.data(), data.size(), MPI_CHAR, dest, data_tag,
                  El::mpi::COMM_WORLD.comm, &requests.back());
      }
    for(const auto index : changed)
      {
        if(owner_root.at(index) != rank || holder.at(index) == rank)
          continue;
        uint64_t size = 0;
        MPI_Recv(&size, 1, MPI_UINT64_T, holder.at(index), size_tag,
                 El::mpi::COMM_WORLD.comm, MPI_STATUS_IGNORE);
        std::string data(size, '\0');
        MPI_Recv(data.data(), size, MPI_CHAR, holder.at(index), data_tag,
                 El::mpi::COMM_WORLD.comm, MPI_STATUS_IGNORE);
        bytes[index] = std::move(data);
      }
    if(!requests.empty())
      MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // Set the local changed blocks, in order (collective in each group)
    for(size_t local = 0; local < block_info.block_indices.size(); ++local)
      {
        const size_t index = block_info.block_indices.at(local);
        if(changed_set.count(index) == 0)
          continue;
        SDP_Block_Data block_data;
        if(is_group_root())
          {
            std::istringstream is(bytes.at(index), std::ios::binary);
            block_data = SDP_Block_Data(is, files.format, local, block_info);
          }
        set_sdp_block_from_root(grid, block_info, block_data, sdp,
                                parameters.solver.sparse_B);
      }
    El::mpi::Barrier(El::mpi::COMM_WORLD);
  }

  void Resident_State::set_objectives_text(const std::string &text)
  {
    std::istringstream is(text);
    read_objectives_from_stream(grid, is, sdp.objective_const,
                                sdp.dual_objective_b);
    files.objectives_text = text;
  }

  // update_pmp: like pmp2sdp --reuseSdp, in memory. "pmp": the new PMP
  // (pmp.json, .nsv, or a list of files). Files listed in the current
  // pmp_info.json keep their blocks; the others are parsed and converted by
  // pmp2sdp's code. The block structure and the normalization must not
  // change. The objective is re-read.
  void Resident_State::cmd_update_pmp(const rapidjson::Value &command,
                                      Json_Writer &w)
  {
    std::vector<fs::path> inputs;
    consistently([&] {
      ASSERT(files.is_directory && !files.block_paths.empty(),
             "update_pmp needs an SDP directory with pmp_info.json");
      for(const auto &s : get_string_list(command, "pmp"))
        inputs.emplace_back(s);
      ASSERT(!inputs.empty(), "update_pmp: no 'pmp' given");
      // Missing files would fail on some ranks only
      if(El::mpi::Rank() == 0)
        for(const auto &file : collect_files_expanding_nsv(inputs))
          ASSERT(fs::exists(file), "update_pmp: no file ", file);
    });

    std::map<fs::path, std::vector<size_t>> blocks_by_path;
    for(size_t index = 0; index < files.block_paths.size(); ++index)
      blocks_by_path[files.block_paths.at(index)].push_back(index);
    const std::function<std::optional<size_t>(const fs::path &)> known
      = [&blocks_by_path](const fs::path &path) -> std::optional<size_t> {
      const auto it = blocks_by_path.find(path);
      if(it == blocks_by_path.end())
        return std::nullopt;
      return it->second.size();
    };

    const auto start = Clock::now();
    Timers timers(env, Verbosity::none);
    std::vector<std::pair<fs::path, size_t>> file_counts;
    // A malformed file fails inside collectives: that aborts, as in pmp2sdp
    auto pmp = read_polynomial_matrix_program(env, inputs, files.max_num_poles,
                                              Verbosity::none, timers, known,
                                              &file_counts);
    const double read_seconds = seconds_since(start);

    // Checks on global data (the same on all ranks), before any change
    std::vector<fs::path> new_paths;
    std::vector<size_t> changed;
    consistently([&] {
      ASSERT_EQUAL(pmp.num_matrices, num_blocks(),
                   "update_pmp: the number of blocks must stay the same");
      for(const auto &[path, count] : file_counts)
        {
          const auto it = blocks_by_path.find(path);
          for(size_t k = 0; k < count; ++k)
            {
              const size_t index = new_paths.size();
              if(it == blocks_by_path.end())
                changed.push_back(index);
              else
                ASSERT(it->second.size() == count
                         && it->second.at(k) == index,
                       "update_pmp: block of ", path,
                       " moved; the block order must stay the same");
              new_paths.push_back(path);
            }
        }
      // The normalization must not change: every block depends on it
      std::string normalization_text;
      if(pmp.normalization.has_value())
        {
          std::ostringstream os;
          write_normalization_json(os, pmp.normalization.value());
          normalization_text = os.str();
        }
      ASSERT(normalization_text == files.normalization_text,
             "update_pmp: the normalization differs from the SDP's");
    });

    const auto convert_start = Clock::now();
    Output_SDP output(pmp, {}, timers);
    std::ostringstream objectives_os;
    write_objectives_json(objectives_os, output.objective_const,
                          output.dual_objective_b);
    consistently([&] {
      ASSERT_EQUAL(static_cast<El::Int>(output.dual_objective_b.size()),
                   sdp.dual_objective_b.Height(),
                   "update_pmp: the dual dimension N must stay the same");
    });
    std::vector<int> holder;
    auto bytes
      = check_new_blocks(output.dual_constraint_groups, changed, holder);

    // Apply
    replace_blocks(std::move(bytes), changed, holder);
    const bool objectives_changed
      = objectives_os.str() != files.objectives_text;
    if(objectives_changed)
      set_objectives_text(objectives_os.str());
    files.block_paths = new_paths;

    w.Key("changed_blocks");
    w.StartArray();
    for(const auto index : changed)
      w.Uint64(index);
    w.EndArray();
    w.Key("objectives_changed");
    w.Bool(objectives_changed);
    w.Key("read_seconds");
    w.Double(read_seconds);
    w.Key("convert_seconds");
    w.Double(seconds_since(convert_start));
  }

  // update_sdp: re-read blocks and objectives from an SDP directory written
  // by pmp2sdp (e.g. with --reuseSdp). "blocks": "changed" (default: the
  // blocks whose source file in pmp_info.json differs), "all", or a list.
  void Resident_State::cmd_update_sdp(const rapidjson::Value &command,
                                      Json_Writer &w)
  {
    fs::path dir;
    consistently([&] {
      dir = get_string(command, "sdp", "");
      ASSERT(!dir.empty(), "update_sdp: no 'sdp' given");
    });
    const auto start = Clock::now();
    Sdp_Files new_files(dir);

    // Validate, the same way on all ranks (the files' text is broadcast)
    std::vector<size_t> changed;
    bool normalization_changed = false;
    consistently([&] {
      ASSERT(new_files.is_directory, "update_sdp: not a directory: ", dir);
      ASSERT_EQUAL(new_files.num_blocks, num_blocks(),
                   "update_sdp: the number of blocks must stay the same");
      if(command.HasMember("blocks") && command["blocks"].IsArray())
        {
          for(const auto &v : command["blocks"].GetArray())
            {
              ASSERT(v.IsUint64() && v.GetUint64() < num_blocks(),
                     "update_sdp: 'blocks' must be block indices below ",
                     num_blocks());
              changed.push_back(v.GetUint64());
            }
          std::sort(changed.begin(), changed.end());
          changed.erase(std::unique(changed.begin(), changed.end()),
                        changed.end());
        }
      else
        {
          const auto mode = get_string(command, "blocks", "changed");
          const bool all = mode == "all";
          ASSERT(all || mode == "changed", "update_sdp: unknown blocks '",
                 mode, "'");
          if(!all)
            ASSERT(!new_files.block_paths.empty()
                     && !files.block_paths.empty(),
                   "update_sdp: blocks 'changed' needs pmp_info.json");
          for(size_t index = 0; index < num_blocks(); ++index)
            if(all
               || new_files.block_paths.at(index)
                    != files.block_paths.at(index))
              changed.push_back(index);
        }
      normalization_changed
        = new_files.normalization_text != files.normalization_text;
      ASSERT(!normalization_changed || changed.size() == num_blocks(),
             "update_sdp: the normalization changed, so all blocks must be "
             "re-read (blocks: \"all\")");
      rapidjson::Document objectives;
      objectives.Parse(new_files.objectives_text.c_str());
      ASSERT(!objectives.HasParseError() && objectives.IsObject()
               && objectives.HasMember("b") && objectives["b"].IsArray()
               && objectives.HasMember("constant"),
             "update_sdp: cannot parse ", dir / "objectives.json");
      ASSERT_EQUAL(static_cast<El::Int>(objectives["b"].Size()),
                   sdp.dual_objective_b.Height(),
                   "update_sdp: the dual dimension N must stay the same");
    });

    // Read and check the changed local blocks on the group roots first
    const auto N = sdp.dual_objective_b.Height();
    std::map<size_t, SDP_Block_Data> new_blocks;
    {
      std::string error;
      if(is_group_root())
        {
          try
            {
              for(size_t local = 0; local < block_info.block_indices.size();
                  ++local)
                {
                  const size_t index = block_info.block_indices.at(local);
                  if(!std::binary_search(changed.begin(), changed.end(),
                                         index))
                    continue;
                  fs::path block_path
                    = dir / ("block_data_" + std::to_string(index) + ".bin");
                  Block_File_Format format = Block_File_Format::bin;
                  if(!fs::exists(block_path))
                    {
                      block_path.replace_extension(".json");
                      format = Block_File_Format::json;
                    }
                  std::ifstream is(block_path, std::ios::binary);
                  ASSERT(is.good(), "Cannot open ", block_path);
                  SDP_Block_Data data(is, format, local, block_info);
                  bool same = data.constraint_matrix.Width() == N
                              && data.constraint_matrix.Height()
                                   == static_cast<El::Int>(
                                     block_info.get_schur_block_size(index))
                              && data.primal_objective_c.Height()
                                   == data.constraint_matrix.Height();
                  for(const size_t parity : {0, 1})
                    same
                      = same
                        && data.bilinear_bases[parity].Height()
                             == static_cast<El::Int>(
                               block_info.get_bilinear_bases_height(index,
                                                                    parity))
                        && data.bilinear_bases[parity].Width()
                             == static_cast<El::Int>(
                               block_info.get_bilinear_bases_width(index,
                                                                   parity));
                  ASSERT(same, "update_sdp: block ", index,
                         " would change its size; the block structure must "
                         "stay the same");
                  new_blocks.emplace(index, std::move(data));
                }
            }
          catch(std::exception &e)
            {
              error = e.what();
            }
        }
      agree_or_throw(error);
    }

    // Apply
    for(size_t local = 0; local < block_info.block_indices.size(); ++local)
      {
        const size_t index = block_info.block_indices.at(local);
        if(!std::binary_search(changed.begin(), changed.end(), index))
          continue;
        SDP_Block_Data empty;
        set_sdp_block_from_root(grid, block_info,
                                is_group_root() ? new_blocks.at(index)
                                                : empty,
                                sdp, parameters.solver.sparse_B);
      }
    // Objectives as SDP's constructor reads them (the same text)
    set_objectives_text(new_files.objectives_text);
    if(normalization_changed)
      {
        normalization.reset();
        if(!new_files.normalization_text.empty())
          {
            rapidjson::Document doc;
            doc.Parse(new_files.normalization_text.c_str());
            std::vector<El::BigFloat> n;
            for(const auto &v : doc["normalization"].GetArray())
              n.emplace_back(v.GetString());
            normalization = n;
          }
        if(El::mpi::Rank() == 0)
          sdp.normalization = normalization;
      }
    El::mpi::Barrier(El::mpi::COMM_WORLD);

    // Source of each block: the new one where it was re-read
    std::vector<fs::path> block_paths(num_blocks());
    for(size_t index = 0; index < num_blocks(); ++index)
      {
        const bool reread
          = std::binary_search(changed.begin(), changed.end(), index);
        if(!reread)
          block_paths.at(index) = files.block_paths.empty()
                                    ? fs::path()
                                    : files.block_paths.at(index);
        else if(!new_files.block_paths.empty())
          block_paths.at(index) = new_files.block_paths.at(index);
        else
          block_paths.at(index) = "unknown:" + dir.string();
      }
    new_files.block_paths = block_paths;
    files = new_files;

    w.Key("changed_blocks");
    w.StartArray();
    for(const auto index : changed)
      w.Uint64(index);
    w.EndArray();
    w.Key("normalization_changed");
    w.Bool(normalization_changed);
    w.Key("read_seconds");
    w.Double(seconds_since(start));
  }

  // load_family: block "block" becomes sum_i w_i (component i). Each
  // component is a PMP file with one matrix of the block's shape (same
  // sample points and bilinear basis), converted with the SDP's
  // normalization. "form_dim": the components are the upper triangle
  // (row-major) of a symmetric form_dim x form_dim matrix, so that
  // set_coefficients can take "lambda" (w_ab = lambda_a lambda_b, doubled
  // off the diagonal). The block keeps its current data until
  // set_coefficients.
  void Resident_State::cmd_load_family(const rapidjson::Value &command,
                                       Json_Writer &w)
  {
    std::string name, normalization_file;
    Family family;
    consistently([&] {
      name = get_string(command, "name", "");
      ASSERT(!name.empty(), "load_family: no 'name' given");
      ASSERT(command.HasMember("block") && command["block"].IsUint64(),
             "load_family: 'block' must be a block index");
      family.block_index = command["block"].GetUint64();
      ASSERT(family.block_index < num_blocks(), "load_family: no block ",
             family.block_index);
      for(const auto &s : get_string_list(command, "components"))
        family.paths.push_back(s);
      family.num_components = family.paths.size();
      ASSERT(family.num_components > 0, "load_family: no components");
      if(command.HasMember("form_dim"))
        {
          ASSERT(command["form_dim"].IsUint64(),
                 "load_family: 'form_dim' must be a positive integer");
          family.form_dim = command["form_dim"].GetUint64();
          ASSERT_EQUAL(family.form_dim * (family.form_dim + 1) / 2,
                       family.num_components,
                       "load_family: form_dim does not match the number of "
                       "components");
        }
      normalization_file = get_string(command, "normalization_file", "");
      ASSERT(normalization.has_value(), "load_family needs a normalization");
    });
    family.owner_root = owner_root.at(family.block_index);
    family.normalization_text = files.normalization_text;
    const auto start = Clock::now();
    std::string error;
    if(El::mpi::Rank() == family.owner_root)
      {
        try
          {
            Timers timers(env, Verbosity::none);
            // The normalization the components are converted with: from the
            // PMP file, parsed as pmp2sdp parses it (exact), or the SDP's
            // normalization.json (a decimal round trip)
            auto conversion_normalization = normalization;
            if(!normalization_file.empty())
              {
                auto parsed = PMP_File_Parse_Result::read(
                  normalization_file, files.max_num_poles, false, true,
                  [](size_t) { return false; });
                ASSERT(parsed.normalization.has_value(), "load_family: ",
                       normalization_file, " has no normalization");
                std::ostringstream os;
                write_normalization_json(os, parsed.normalization.value());
                ASSERT(os.str() == files.normalization_text, "load_family: ",
                       normalization_file,
                       " has another normalization than the SDP");
                conversion_normalization = parsed.normalization;
              }
            const size_t N1 = conversion_normalization->size();
            for(const auto &path : family.paths)
              {
                auto parsed = PMP_File_Parse_Result::read(
                  path, files.max_num_poles, false, false,
                  [](size_t) { return true; });
                ASSERT_EQUAL(parsed.num_matrices, 1, path);
                auto &matrix = parsed.parsed_matrices.at(0);
                ASSERT_EQUAL(matrix.polynomials(0, 0).size(), N1, path);
                family.first_scaling.push_back(matrix.sample_scalings.at(0));
                std::vector<Polynomial_Vector_Matrix> matrices;
                matrices.emplace_back(std::move(matrix));
                Polynomial_Matrix_Program pmp(
                  std::vector<El::BigFloat>(N1, El::BigFloat(0)),
                  conversion_normalization, num_blocks(), std::move(matrices),
                  {family.block_index}, {path});
                Output_SDP output(pmp, {}, timers);
                const auto &group = output.dual_constraint_groups.at(0);
                ASSERT(group.dim
                           == block_info.dimensions.at(family.block_index)
                         && group.num_points
                              == block_info.num_points.at(family.block_index),
                       "load_family: ", path,
                       " does not have the shape of block ",
                       family.block_index);
                // As if read from a block file
                std::ostringstream os(std::ios::binary);
                write_block_data(os, group, files.format);
                std::istringstream is(os.str(), std::ios::binary);
                auto data = parse_block_data(is, files.format);
                if(family.c.empty())
                  {
                    family.bilinear_bases[0] = data.bilinear_bases_even;
                    family.bilinear_bases[1] = data.bilinear_bases_odd;
                  }
                else
                  {
                    for(const size_t parity : {0, 1})
                      {
                        const auto &a = family.bilinear_bases[parity];
                        const auto &b = parity == 0 ? data.bilinear_bases_even
                                                    : data.bilinear_bases_odd;
                        bool equal = a.Height() == b.Height()
                                     && a.Width() == b.Width();
                        for(int i = 0; equal && i < a.Height(); ++i)
                          for(int j = 0; equal && j < a.Width(); ++j)
                            equal = a(i, j) == b(i, j);
                        ASSERT(equal, "load_family: ", path,
                               " has a different bilinear basis");
                      }
                  }
                El::Matrix<El::BigFloat> c(data.c.size(), 1);
                for(size_t i = 0; i < data.c.size(); ++i)
                  c(i, 0) = data.c.at(i);
                family.c.push_back(std::move(c));
                family.B.push_back(std::move(data.B));
              }
          }
        catch(std::exception &e)
          {
            error = e.what();
          }
      }
    agree_or_throw(error);
    families[name] = std::move(family);
    w.Key("name");
    w.String(name.c_str());
    w.Key("block");
    w.Uint64(families[name].block_index);
    w.Key("components");
    w.Uint64(families[name].num_components);
    w.Key("read_seconds");
    w.Double(seconds_since(start));
  }

  // set_coefficients: set the family's block to sum_i w_i (c_i, B_i).
  // "coefficients": w_i (decimal strings or p/q), or "lambda": the vector of
  // a symmetric family (w_ab = lambda_a lambda_b, times 2 for a < b).
  void Resident_State::cmd_set_coefficients(const rapidjson::Value &command,
                                            Json_Writer &w)
  {
    std::string name;
    std::vector<El::BigFloat> coefficients;
    consistently([&] {
      name = get_string(command, "family", "");
      const auto it = families.find(name);
      ASSERT(it != families.end(), "set_coefficients: no family '", name,
             "'");
      const auto &f = it->second;
      ASSERT(f.normalization_text == files.normalization_text,
             "set_coefficients: family '", name,
             "' was loaded with another normalization; load it again");
      if(command.HasMember("lambda"))
        {
          ASSERT(f.form_dim > 0,
                 "set_coefficients: 'lambda' needs a family with form_dim");
          std::vector<El::BigFloat> lambda;
          for(const auto &s : get_string_list(command, "lambda"))
            lambda.push_back(parse_number(s));
          ASSERT_EQUAL(lambda.size(), f.form_dim,
                       "set_coefficients: wrong number of lambdas");
          for(size_t a = 0; a < f.form_dim; ++a)
            for(size_t b = a; b < f.form_dim; ++b)
              coefficients.push_back(a == b ? lambda.at(a) * lambda.at(b)
                                            : 2 * lambda.at(a) * lambda.at(b));
        }
      else
        {
          for(const auto &s : get_string_list(command, "coefficients"))
            coefficients.push_back(parse_number(s));
        }
      ASSERT_EQUAL(coefficients.size(), f.num_components,
                   "set_coefficients: wrong number of coefficients");
    });
    const auto &family = families.at(name);

    const auto start = Clock::now();
    const auto local = local_index(family.block_index);
    if(local.has_value())
      {
        SDP_Block_Data block_data;
        if(is_group_root())
          {
            block_data.block_index_local = local.value();
            auto c = family.c.at(0);
            auto B = family.B.at(0);
            El::Scale(coefficients.at(0), c);
            El::Scale(coefficients.at(0), B);
            for(size_t i = 1; i < family.num_components; ++i)
              {
                El::Axpy(coefficients.at(i), family.c.at(i), c);
                El::Axpy(coefficients.at(i), family.B.at(i), B);
              }
            block_data.primal_objective_c = std::move(c);
            block_data.constraint_matrix = std::move(B);
            block_data.bilinear_bases = family.bilinear_bases;
            for(const size_t parity : {0, 1})
              {
                auto &bases_block = block_data.bases_blocks[parity];
                bases_block.Resize(
                  block_info.get_psd_matrix_block_size(family.block_index,
                                                       parity),
                  block_info.get_bilinear_pairing_block_size(
                    family.block_index, parity));
                set_bilinear_bases_block_local(
                  block_data.bilinear_bases[parity], bases_block);
              }
          }
        set_sdp_block_from_root(grid, block_info, block_data, sdp,
                                parameters.solver.sparse_B);
      }
    El::mpi::Barrier(El::mpi::COMM_WORLD);
    // The block no longer comes from a file
    if(!files.block_paths.empty())
      files.block_paths.at(family.block_index)
        = "family:" + name + ":" + std::to_string(num_solves);
    w.Key("family");
    w.String(name.c_str());
    w.Key("set_seconds");
    w.Double(seconds_since(start));
  }

  // negate_normalization: n -> -n. With the same normalization index, this
  // negates c in every block and the objective constant, and leaves B and b
  // unchanged, exactly (sign-magnitude arithmetic). For example the upper vs
  // lower bound of an OPE coefficient at a fixed point.
  void Resident_State::cmd_negate_normalization(Json_Writer &w)
  {
    consistently([&] {
      ASSERT(normalization.has_value(),
             "negate_normalization: no normalization");
    });
    for(auto &block : sdp.primal_objective_c.blocks)
      {
        auto &local = block.Matrix();
        for(int i = 0; i < local.Height(); ++i)
          for(int j = 0; j < local.Width(); ++j)
            local(i, j) = -local(i, j);
      }
    sdp.objective_const = -sdp.objective_const;
    for(auto &x : normalization.value())
      x = -x;
    if(El::mpi::Rank() == 0)
      sdp.normalization = normalization;
    std::ostringstream os;
    write_normalization_json(os, normalization.value());
    // Families with the old normalization follow: their c_i change sign too
    for(auto &[family_name, family] : families)
      {
        if(family.normalization_text != files.normalization_text)
          continue;
        for(auto &c : family.c) // on the owner root only
          for(int i = 0; i < c.Height(); ++i)
            for(int j = 0; j < c.Width(); ++j)
              c(i, j) = -c(i, j);
        family.normalization_text = os.str();
      }
    files.normalization_text = os.str();
    // The objective constant changed sign: the next update_pmp sets the
    // objective from its files in any case
    files.objectives_text.clear();
    El::mpi::Barrier(El::mpi::COMM_WORLD);
    w.Key("objective_const");
    write_bigfloat(w, sdp.objective_const);
  }

  // save_iterate: "name", in memory; or with "dir", as a checkpoint there
  void Resident_State::cmd_save_iterate(const rapidjson::Value &command,
                                        Json_Writer &w)
  {
    std::string name, dir;
    consistently([&] {
      name = get_string(command, "name", "");
      dir = get_string(command, "dir", "");
      ASSERT(!name.empty() || !dir.empty(),
             "save_iterate: give 'name' or 'dir'");
    });
    if(!dir.empty())
      {
        // A separate checkpoint: keep the counters of the main one
        checkpoint_elsewhere(dir);
        w.Key("dir");
        w.String(dir.c_str());
      }
    else
      {
        iterates.erase(name);
        iterates.emplace(name,
                         Iterate{solver.x, solver.y, solver.X, solver.Y});
        w.Key("name");
        w.String(name.c_str());
      }
  }

  // restore_iterate: "name" (in memory), or "dir" (a checkpoint)
  void Resident_State::cmd_restore_iterate(const rapidjson::Value &command,
                                           Json_Writer &w)
  {
    std::string name, dir;
    consistently([&] {
      name = get_string(command, "name", "");
      dir = get_string(command, "dir", "");
      ASSERT(!name.empty() || !dir.empty(),
             "restore_iterate: give 'name' or 'dir'");
      if(dir.empty())
        ASSERT(iterates.count(name) != 0, "restore_iterate: no iterate '",
               name, "'");
    });
    if(!dir.empty())
      {
        check_checkpoint_dir(dir);
        const auto current = solver.current_generation;
        const auto backup = solver.backup_generation;
        solver.load_checkpoint(dir, block_info, parameters.verbosity, true);
        solver.current_generation = current;
        solver.backup_generation = backup;
      }
    else
      {
        const auto &it = iterates.at(name);
        solver.x = it.x;
        solver.X = it.X;
        solver.y = it.y;
        solver.Y = it.Y;
      }
    dirty = true;
    w.Key("restored");
    w.String(dir.empty() ? name.c_str() : dir.c_str());
  }

  void Resident_State::cmd_status(Json_Writer &w)
  {
    w.Key("solves");
    w.Uint64(num_solves);
    w.Key("sdpDir");
    w.String(files.dir.string().c_str());
    w.Key("iterates");
    w.StartArray();
    for(const auto &[name, iterate] : iterates)
      w.String(name.c_str());
    w.EndArray();
    w.Key("families");
    w.StartArray();
    for(const auto &[name, family] : families)
      w.String(name.c_str());
    w.EndArray();
  }
}

void run_resident(const Environment &env, const SDPB_Parameters &parameters,
                  const std::vector<std::string> &startup_args,
                  const mp_bitcnt_t startup_precision,
                  const Block_Info &block_info)
{
  const auto start = Clock::now();
  Control_Channel channel(parameters.control_dir,
                          parameters.control_idle_timeout);
  Resident_State state(env, parameters, startup_args, startup_precision,
                       block_info);
  if(El::mpi::Rank() == 0)
    {
      channel.write_file("ready.json", state.ready_json(seconds_since(start)));
      if(parameters.verbosity >= Verbosity::regular)
        El::Output(boost::posix_time::second_clock::local_time(),
                   " Resident mode: waiting for commands in ",
                   parameters.control_dir);
    }

  bool keep_running = true;
  while(keep_running)
    {
      const auto text = channel.next_command(env);
      if(!text.has_value())
        break;
      rapidjson::Document document;
      document.Parse(text->c_str());
      rapidjson::StringBuffer sb;
      Json_Writer w(sb);
      w.StartObject();
      w.Key("seq");
      w.Uint64(channel.seq());
      if(document.HasParseError())
        {
          w.Key("error");
          w.String("cannot parse the command as JSON");
        }
      else
        {
          // One command, or {"commands": [...]} executed in order until the
          // first error
          std::vector<const rapidjson::Value *> commands;
          const bool batch
            = document.IsObject() && document.HasMember("commands");
          if(batch && document["commands"].IsArray())
            for(const auto &c : document["commands"].GetArray())
              commands.push_back(&c);
          else if(!batch)
            commands.push_back(&document);
          w.Key("results");
          w.StartArray();
          if(batch && !document["commands"].IsArray())
            {
              w.StartObject();
              w.Key("error");
              w.String("'commands' must be an array");
              w.EndObject();
            }
          for(const auto *command : commands)
            {
              rapidjson::StringBuffer csb;
              Json_Writer cw(csb);
              std::string error;
              try
                {
                  keep_running = state.execute(*command, cw);
                }
              catch(Command_Error &e)
                {
                  // The same on all ranks, before any change: carry on
                  error = e.what();
                }
              catch(std::exception &e)
                {
                  // Possibly on some ranks only, inside a collective
                  El::ReportException(e);
                  El::mpi::Abort(El::mpi::COMM_WORLD, 1);
                }
              if(error.empty())
                {
                  w.RawValue(csb.GetString(), csb.GetSize(),
                             rapidjson::kObjectType);
                }
              else
                {
                  w.StartObject();
                  w.Key("error");
                  w.String(error.c_str());
                  w.EndObject();
                  if(El::mpi::Rank() == 0)
                    El::Output("Resident mode: command error: ", error);
                  break;
                }
              if(!keep_running)
                break;
            }
          w.EndArray();
        }
      w.EndObject();
      if(El::mpi::Rank() == 0)
        channel.write_result(sb.GetString());
    }

  // Stopped without 'quit' (SIGTERM, stop file, idle timeout): keep the
  // iterate if it is newer than the last checkpoint
  auto reason = channel.stop_reason();
  if(state.sigterm)
    reason = Control_Channel::Stop_Reason::sigterm;
  if(reason != Control_Channel::Stop_Reason::none && state.dirty)
    state.checkpoint(parameters.solver.checkpoint_out,
                     to_property_tree(parameters));
  if(El::mpi::Rank() == 0)
    {
      rapidjson::StringBuffer sb;
      Json_Writer w(sb);
      w.StartObject();
      w.Key("status");
      w.String("exit");
      w.Key("reason");
      const auto reason_string
        = reason == Control_Channel::Stop_Reason::none ? std::string("quit")
                                                       : to_string(reason);
      w.String(reason_string.c_str());
      w.Key("last_seq");
      w.Uint64(channel.seq());
      w.Key("pid");
      w.Int64(getpid());
      w.Key("seconds");
      w.Double(seconds_since(start));
      w.EndObject();
      channel.write_file("exit.json", sb.GetString());
    }
  // As sdpb does after a SIGTERM: exit code SIGTERM, so that job scripts
  // that requeue on it see it
  if(reason == Control_Channel::Stop_Reason::sigterm)
    {
      if(El::mpi::Rank() == 0)
        El::Output("Received SIGTERM, exiting gracefully...");
      MPI_Finalize();
      std::exit(SIGTERM);
    }
}

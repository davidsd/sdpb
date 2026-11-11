#include "Reuse_SDP.hxx"

#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <rapidjson/document.h>

#include <fstream>
#include <set>
#include <sstream>

namespace fs = std::filesystem;

// Boost's default command line style also accepts unambiguous prefixes
// of long names ("--maxNum 2"). A command line that abbreviates one of the
// options read here is refused. pmp2sdp's own long options are exact
// names, not abbreviations ("--output" is not "--outputFormat").
static void
refuse_abbreviation(const std::string &t, const std::string &name)
{
  if(t.rfind("--", 0) != 0)
    return;
  const std::string key = t.substr(0, t.find('='));
  static const std::set<std::string> long_options{
    "--help",         "--input",    "--output",    "--precision",
    "--maxNumPoles",  "--outputFormat",            "--zip",
    "--reuseSdp",     "--verbosity", "--version"};
  if(long_options.count(key) != 0)
    return;
  if(key.size() > 2 && key.size() < name.size()
     && name.rfind(key, 0) == 0)
    RUNTIME_ERROR("cannot tell the SDP's ", name,
                  " from the abbreviated option ", t,
                  " in the command its control.json records; write it with "
                  "full option names.");
}

std::optional<std::string>
find_pmp2sdp_option(const std::vector<std::string> &tokens,
                    const std::string &name, const std::string &short_name)
{
  for(size_t i = 0; i < tokens.size(); ++i)
    {
      const auto &t = tokens.at(i);
      refuse_abbreviation(t, name);
      if(t == name || t == short_name)
        {
          ASSERT(i + 1 < tokens.size(), "Missing value for ", t);
          return tokens.at(i + 1);
        }
      if(t.rfind(name + "=", 0) == 0)
        return t.substr(name.size() + 1);
      if(t.size() > short_name.size() && t.rfind(short_name, 0) == 0)
        return t.substr(short_name.size());
    }
  return std::nullopt;
}

namespace
{
  std::string read_file(const fs::path &path)
  {
    std::ifstream in(path, std::ios::binary);
    ASSERT(in.good(), "Failed to open ", path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  void broadcast_string(std::string &s)
  {
    size_t size = s.size();
    El::mpi::Broadcast(size, 0, El::mpi::COMM_WORLD);
    s.resize(size);
    if(size > 0)
      MPI_Bcast(s.data(), size, MPI_CHAR, 0, El::mpi::COMM_WORLD.comm);
  }

  bool has_flag(const std::vector<std::string> &tokens,
                const std::string &name, const std::string &short_name)
  {
    for(const auto &t : tokens)
      {
        refuse_abbreviation(t, name);
        if(t == name || t == short_name || t.rfind(name + "=", 0) == 0)
          return true;
      }
    return false;
  }

  // Check that the old SDP was written with the same conversion settings.
  void check_old_command(const fs::path &old_sdp, const std::string &control,
                         const int precision, const int64_t max_num_poles,
                         const Block_File_Format format)
  {
    rapidjson::Document document;
    document.Parse(control.c_str());
    ASSERT(!document.HasParseError() && document.HasMember("command"),
           "Cannot read command from ", old_sdp / "control.json");
    std::vector<std::string> tokens;
    {
      std::istringstream iss(document["command"].GetString());
      std::string token;
      while(iss >> token)
        tokens.push_back(token);
    }
    const auto old_precision
      = find_pmp2sdp_option(tokens, "--precision", "-p");
    ASSERT(old_precision.has_value(),
           "Cannot determine --precision of the old SDP from ",
           old_sdp / "control.json");
    ASSERT_EQUAL(std::stoi(old_precision.value()), precision,
                 "--reuseSdp requires the old SDP to have the same precision");
    const auto old_poles
      = find_pmp2sdp_option(tokens, "--maxNumPoles", "-n");
    ASSERT_EQUAL(old_poles.has_value() ? std::stoll(old_poles.value()) : -1,
                 max_num_poles,
                 "--reuseSdp requires the old SDP to have the same "
                 "--maxNumPoles");
    const auto old_format
      = find_pmp2sdp_option(tokens, "--outputFormat", "-f");
    const std::string format_name
      = format == Block_File_Format::bin ? "bin" : "json";
    ASSERT_EQUAL(old_format.value_or("bin"), format_name,
                 "--reuseSdp requires the old SDP to have the same "
                 "--outputFormat");
    ASSERT(!has_flag(tokens, "--zip", "-z"),
           "--reuseSdp does not support a zipped old SDP");
  }

  // Split a JSON array of objects, as written by write_pmp_info_json(),
  // into the raw text of its elements.
  std::vector<std::string> split_json_array(const std::string &text)
  {
    std::vector<std::string> result;
    size_t pos = text.find('[');
    ASSERT(pos != std::string::npos, "pmp_info.json: '[' not found");
    ++pos;
    int depth = 0;
    bool in_string = false;
    size_t begin = 0;
    for(; pos < text.size(); ++pos)
      {
        const char c = text[pos];
        if(in_string)
          {
            if(c == '\\')
              ++pos;
            else if(c == '"')
              in_string = false;
            continue;
          }
        if(c == '"')
          in_string = true;
        else if(c == '{' || c == '[')
          {
            if(depth == 0)
              {
                ASSERT(c == '{', "pmp_info.json: expected an object");
                begin = pos;
              }
            ++depth;
          }
        else if(c == '}' || c == ']')
          {
            if(depth == 0)
              {
                ASSERT(c == ']', "pmp_info.json: unbalanced braces");
                return result;
              }
            --depth;
            if(depth == 0)
              result.push_back(text.substr(begin, pos + 1 - begin));
          }
      }
    RUNTIME_ERROR("pmp_info.json: unterminated array");
  }

  const std::string index_prefix = "{\"index\":";
}

Reused_SDP::Reused_SDP(const fs::path &old_sdp, const int precision,
                       const int64_t max_num_poles,
                       const Block_File_Format format)
    : old_sdp(old_sdp), format(format)
{
  std::string control, pmp_info;
  if(El::mpi::Rank() == 0)
    {
      ASSERT(fs::is_directory(old_sdp),
             "--reuseSdp: not a directory: ", old_sdp);
      control = read_file(old_sdp / "control.json");
      pmp_info = read_file(old_sdp / "pmp_info.json");
    }
  broadcast_string(control);
  broadcast_string(pmp_info);

  check_old_command(old_sdp, control, precision, max_num_poles, format);

  rapidjson::Document control_doc;
  control_doc.Parse(control.c_str());
  const size_t old_num_blocks = control_doc["num_blocks"].GetInt();

  raw_pmp_info_entries = split_json_array(pmp_info);
  ASSERT_EQUAL(raw_pmp_info_entries.size(), old_num_blocks,
               "Number of entries in ", old_sdp / "pmp_info.json");
  for(size_t i = 0; i < raw_pmp_info_entries.size(); ++i)
    {
      const auto &entry = raw_pmp_info_entries.at(i);
      ASSERT(entry.rfind(index_prefix, 0) == 0,
             "pmp_info.json entry does not start with ", index_prefix);
      rapidjson::Document doc;
      doc.Parse(entry.c_str());
      ASSERT(!doc.HasParseError(), "Cannot parse pmp_info.json entry ", i);
      ASSERT_EQUAL(static_cast<size_t>(doc["index"].GetInt()), i,
                   "pmp_info.json entries must be ordered by index");
      old_blocks_by_path[fs::path(doc["path"].GetString())].push_back(i);
    }
  // Only rank 0 writes pmp_info.json
  if(El::mpi::Rank() != 0)
    {
      for(auto &entry : raw_pmp_info_entries)
        std::string().swap(entry);
    }
}

std::optional<size_t>
Reused_SDP::num_matrices(const fs::path &input_path) const
{
  const auto it = old_blocks_by_path.find(input_path);
  if(it == old_blocks_by_path.end())
    return std::nullopt;
  return it->second.size();
}

void Reused_SDP::set_new_block_mapping(
  const std::vector<std::pair<fs::path, size_t>> &files)
{
  new_to_old.clear();
  std::set<fs::path> reused_paths;
  for(const auto &[path, count] : files)
    {
      const auto it = old_blocks_by_path.find(path);
      if(it != old_blocks_by_path.end())
        {
          ASSERT(reused_paths.insert(path).second,
                 "--reuseSdp: input file listed more than once: ", path);
          ASSERT_EQUAL(it->second.size(), count, DEBUG_STRING(path));
        }
      for(size_t j = 0; j < count; ++j)
        {
          if(it == old_blocks_by_path.end())
            new_to_old.emplace_back(std::nullopt);
          else
            new_to_old.emplace_back(it->second.at(j));
        }
    }
}

size_t Reused_SDP::num_reused_blocks() const
{
  size_t result = 0;
  for(const auto &old : new_to_old)
    if(old.has_value())
      ++result;
  return result;
}

std::string Reused_SDP::renumbered_pmp_info_entry(const size_t old_index,
                                                  const size_t new_index) const
{
  const auto &entry = raw_pmp_info_entries.at(old_index);
  const size_t end = entry.find(',', index_prefix.size());
  ASSERT(end != std::string::npos);
  return index_prefix + std::to_string(new_index) + entry.substr(end);
}

fs::path Reused_SDP::normalization_path() const
{ return old_sdp / "normalization.json"; }

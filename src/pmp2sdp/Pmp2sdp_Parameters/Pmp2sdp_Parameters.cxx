#include "Pmp2sdp_Parameters.hxx"

#include "sdpb_util/assert.hxx"

#include <El.hpp>

namespace po = boost::program_options;
namespace fs = std::filesystem;

Pmp2sdp_Parameters::Pmp2sdp_Parameters(int argc, char **argv)
{
  for(int arg(0); arg != argc; ++arg)
    {
      command_arguments.emplace_back(argv[arg]);
    }

  po::options_description options("Basic options");
  options.add_options()("help,h", "Show this helpful message.");
  options.add_options()(
    "input,i", po::value<fs::path>(&input_file)->required(),
    "Mathematica, JSON, XML or NSV file with SDP definition");
  options.add_options()("output,o",
                        po::value<fs::path>(&output_path)->required(),
                        "Directory to place output");
  options.add_options()(
    "precision,p", po::value<int>(&precision)->required(),
    "The precision, in the number of bits, for numbers in the "
    "computation. ");
  options.add_options()(
    "maxNumPoles,n", po::value<int64_t>(&max_num_poles)->default_value(-1),
    "Maximum number of poles in reducedPrefactor (-1 means unlimited).\n"
    "pmp2sdp will keep up to maxNumPoles rightmost poles in reducedPrefactor, "
    "which is used for interpolating the PMP.\n"
    "See https://arxiv.org/abs/2509.14307, Section 4.3 for details.");
  options.add_options()(
    "encodingPoints",
    po::value<std::string>(&encoding_points)->default_value("prefactorMatched"),
    "Two-set encoding. The polynomials stay defined at the approximation "
    "points (the sample points: arXiv:2509.14307's nodes, or the input's), "
    "and the SDP is written at a second set of as many encoding points, "
    "from the values of the interpolant at the approximation points, with "
    "the same scalings w_j for B, c and the bilinear bases (made by the "
    "three-term recurrence). The polynomial program, so the feasible set "
    "and the optimum, are unchanged. Rules, with a = encodingPointsScale * "
    "median |kept pole| (10 * encodingPointsScale without poles) or a from "
    "--encodingScaleFile: 'functionalMatched': the quadrature rule exact "
    "to degree d for dx / (sqrt(x) g(x)), g = prod (x - p_i) "
    "(x + a)^(d + 1 - #poles) (points at N_g(y_k) = k - 1/2, scalings "
    "1/P_g; see docs/Two_Set_Encoding.md); 'chebyshev': "
    "y_k = a tan^2((2k-1) pi / (4(d+1))), scalings 1/(prod (y_k - p_i) "
    "(y_k + a)^(d - #poles)); 'file': points and scalings from "
    "--encodingPointsFile. 'prefactorMatched' (default): the SDP is written "
    "at the approximation points, as without this option.");
  options.add_options()(
    "encodingPointsScale",
    po::value<double>(&encoding_points_scale)->default_value(1),
    "With --encodingPoints functionalMatched or chebyshev: a = this factor "
    "times the median of |p_i| over the block's kept poles (times 10 for a "
    "block without poles).");
  options.add_options()(
    "encodingScaleFile", po::value<fs::path>(&encoding_scale_file),
    "With --encodingPoints functionalMatched or chebyshev: a JSON array of "
    "the scale a of every block, by block index (replaces "
    "--encodingPointsScale).");
  options.add_options()(
    "encodingPointsFile", po::value<fs::path>(&encoding_points_file),
    "With --encodingPoints file: a JSON file with the encoding points and "
    "scalings of every block, by block index: either a PMP-like object "
    "{\"PositiveMatrixWithPrefactorArray\": [{\"samplePoints\": [...], "
    "\"sampleScalings\": [...]}, ...]} or the array alone. Other fields "
    "are ignored.");
  options.add_options()(
    "encodingExtraPrecision",
    po::value<int64_t>(&encoding_extra_precision)->default_value(384),
    "With --encodingPoints: the conversion runs at precision + this many "
    "bits, and the SDP is written at --precision. The step to the encoding "
    "points cancels digits (printed); pmp2sdp stops if fewer than 32 bits "
    "are left over. (Measured on 3d Ising bootstrap programs: 45 to 71 "
    "digits.)");
  options.add_options()(
    "outputFormat,f",
    po::value<Block_File_Format>(&output_format)
      ->default_value(Block_File_Format::bin),
    "Output format for SDP blocks. Could be either 'bin' or 'json'");
  options.add_options()(
    "zip,z", po::bool_switch(&zip),
    "Store output to zip file instead of plain directory.");
  options.add_options()(
    "verbosity,v",
    po::value<Verbosity>(&verbosity)->default_value(Verbosity::regular),
    "Verbosity.  0 -> no output, 1 -> regular output, 2 -> debug output, 3 -> "
    "trace output");

  // (partial) backward compatibility with pvm2sdp
  po::positional_options_description positional;
  positional.add("precision", 1);
  positional.add("input", 1);
  positional.add("output", 1);

  try
    {
      po::variables_map variables_map;

      po::store(po::command_line_parser(argc, argv)
                  .options(options)
                  .positional(positional)
                  .run(),
                variables_map);

      if(variables_map.count("help") != 0)
        {
          if(El::mpi::Rank() == 0)
            El::Output(options);
          return;
        }

      po::notify(variables_map);

      ASSERT(fs::exists(input_file),
             "Input file does not exist: ", input_file);
      ASSERT(!fs::is_directory(input_file) && input_file != ".",
             "Input file is a directory, not a file:", input_file);
      ASSERT(encoding_points == "prefactorMatched"
               || encoding_points == "functionalMatched"
               || encoding_points == "chebyshev" || encoding_points == "file",
             "--encodingPoints must be prefactorMatched, functionalMatched, "
             "chebyshev or file: ",
             encoding_points);
      if(two_set_encoding())
        {
          ASSERT(encoding_extra_precision >= 64,
                 "--encodingExtraPrecision must be at least 64 bits");
        }
      ASSERT(encoding_points_scale > 0,
             "--encodingPointsScale must be positive");
      ASSERT((encoding_points == "file") == !encoding_points_file.empty(),
             "--encodingPointsFile goes with --encodingPoints file, and "
             "only with it");
      ASSERT(encoding_scale_file.empty()
               || encoding_points == "functionalMatched"
               || encoding_points == "chebyshev",
             "--encodingScaleFile needs --encodingPoints functionalMatched "
             "or chebyshev");
    }
  catch(po::error &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
bool Pmp2sdp_Parameters::is_valid() const
{ return !input_file.empty(); }
boost::property_tree::ptree to_property_tree(const Pmp2sdp_Parameters &p)
{
  boost::property_tree::ptree result;

  result.put("input", p.input_file.string());
  result.put("output", p.output_path.string());
  result.put("precision", p.precision);
  result.put("maxNumPoles", p.max_num_poles);
  if(p.two_set_encoding())
    {
      result.put("encodingPoints", p.encoding_points);
      result.put("encodingPointsScale", p.encoding_points_scale);
      result.put("encodingScaleFile", p.encoding_scale_file.string());
      result.put("encodingPointsFile", p.encoding_points_file.string());
      result.put("encodingExtraPrecision", p.encoding_extra_precision);
    }
  result.put("outputFormat", p.output_format);
  result.put("zip", p.zip);
  result.put("verbosity", static_cast<int>(p.verbosity));

  return result;
}
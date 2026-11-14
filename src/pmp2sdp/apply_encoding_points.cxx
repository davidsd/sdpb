// pmp2sdp --encodingPoints: set every block's encoding points and scalings
// (see pmp/convert/encoding_points.cxx and
// Dual_Constraint_Group/two_set_encoding.cxx).

#include "apply_encoding_points.hxx"

#include "sdpb_util/assert.hxx"

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

#include <fstream>

namespace fs = std::filesystem;

namespace
{
  rapidjson::Document read_json(const fs::path &path)
  {
    std::ifstream input(path);
    ASSERT(input.good(), "Cannot open ", path);
    rapidjson::IStreamWrapper wrapper(input);
    rapidjson::Document document;
    document.ParseStream<rapidjson::kParseNumbersAsStringsFlag>(wrapper);
    ASSERT(!document.HasParseError(), "JSON parse error in ", path,
           " at offset ", document.GetErrorOffset());
    return document;
  }

  El::BigFloat to_BigFloat(const rapidjson::Value &value,
                           const fs::path &path)
  {
    ASSERT(value.IsString(), "expected a number or a string in ", path);
    return El::BigFloat(std::string(value.GetString(),
                                    value.GetStringLength()));
  }

  std::vector<El::BigFloat> to_BigFloat_Vector(const rapidjson::Value &value,
                                               const fs::path &path)
  {
    ASSERT(value.IsArray(), "expected an array in ", path);
    std::vector<El::BigFloat> result;
    for(const auto &element : value.GetArray())
      result.push_back(to_BigFloat(element, path));
    return result;
  }

  const rapidjson::Value &
  block_array(const rapidjson::Document &document, const fs::path &path)
  {
    if(document.IsObject())
      {
        ASSERT(document.HasMember("PositiveMatrixWithPrefactorArray"),
               "--encodingPointsFile: no PositiveMatrixWithPrefactorArray in ",
               path);
        const auto &array = document["PositiveMatrixWithPrefactorArray"];
        ASSERT(array.IsArray(), path);
        return array;
      }
    ASSERT(document.IsArray(),
           "--encodingPointsFile: expected an object or an array in ", path);
    return document;
  }
}

Encoding_Points_Report
apply_encoding_points(Polynomial_Matrix_Program &pmp,
                      const Pmp2sdp_Parameters &parameters)
{
  Encoding_Points_Report report;
  if(!parameters.two_set_encoding())
    return report;

  std::optional<rapidjson::Document> scale_document, points_document;
  if(!parameters.encoding_scale_file.empty())
    {
      scale_document = read_json(parameters.encoding_scale_file);
      ASSERT(scale_document->IsArray(),
             "--encodingScaleFile: expected an array of numbers in ",
             parameters.encoding_scale_file);
    }
  if(parameters.encoding_points == "file")
    points_document = read_json(parameters.encoding_points_file);

  const El::BigFloat scale_factor(parameters.encoding_points_scale);
  for(size_t i = 0; i < pmp.matrices.size(); ++i)
    {
      auto &m = pmp.matrices.at(i);
      const size_t block_index = pmp.matrix_index_local_to_global.at(i);
      const size_t num_points = m.sample_points.size();
      // A degree-0 block (one sample point, e.g. an isolated constant
      // constraint) has nothing to re-encode: its single approximation point
      // is its encoding.  Leave it on the default path (the reach ratio below
      // would also divide by its point x = 0); main() then brings it to the
      // re-encoded blocks' scale (scale_degree_0_blocks).
      if(num_points < 2)
        continue;
      std::vector<El::BigFloat> points, scalings;
      if(points_document.has_value())
        {
          const auto &path = parameters.encoding_points_file;
          const auto &array = block_array(points_document.value(), path);
          ASSERT(block_index < array.Size(), "--encodingPointsFile ", path,
                 " has ", array.Size(), " blocks, block ", block_index,
                 " is missing");
          const auto &block = array[block_index];
          ASSERT(block.IsObject() && block.HasMember("samplePoints")
                   && block.HasMember("sampleScalings"),
                 "--encodingPointsFile: block ", block_index,
                 " needs samplePoints and sampleScalings");
          points = to_BigFloat_Vector(block["samplePoints"], path);
          scalings = to_BigFloat_Vector(block["sampleScalings"], path);
        }
      else
        {
          El::BigFloat a;
          if(scale_document.has_value())
            {
              const auto &path = parameters.encoding_scale_file;
              ASSERT(block_index < scale_document->Size(),
                     "--encodingScaleFile ", path, " has ",
                     scale_document->Size(), " entries, block ", block_index,
                     " is missing");
              a = to_BigFloat((*scale_document)[block_index], path);
              ASSERT(a > 0, "--encodingScaleFile: a must be positive",
                     DEBUG_STRING(block_index));
            }
          else
            a = encoding_points_scale(m.reduced_prefactor, scale_factor);

          if(parameters.encoding_points == "functionalMatched")
            functional_matched_encoding_points(m.reduced_prefactor,
                                               num_points, a, points,
                                               scalings);
          else if(parameters.encoding_points == "chebyshev")
            chebyshev_encoding_points(m.reduced_prefactor, num_points, a,
                                      points, scalings);
          else
            RUNTIME_ERROR("unknown --encodingPoints ",
                          parameters.encoding_points);
        }
      // Largest encoding point over the largest approximation point
      const El::BigFloat ratio
        = *std::max_element(points.begin(), points.end())
          / *std::max_element(m.sample_points.begin(), m.sample_points.end());
      report.max_reach = std::max(report.max_reach, static_cast<double>(ratio));
      set_encoding_points(m, std::move(points), std::move(scalings));
    }
  return report;
}

void truncate_precision(El::BigFloat &x, const mp_bitcnt_t precision)
{
  x.SetPrecision(precision);
}

void truncate_precision(std::vector<El::BigFloat> &v,
                        const mp_bitcnt_t precision)
{
  for(auto &x : v)
    x.SetPrecision(precision);
}

void truncate_precision(El::Matrix<El::BigFloat> &m,
                        const mp_bitcnt_t precision)
{
  for(El::Int j = 0; j < m.Width(); ++j)
    for(El::Int i = 0; i < m.Height(); ++i)
      m.Ref(i, j).SetPrecision(precision);
}

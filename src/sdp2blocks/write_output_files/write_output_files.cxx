#include "../Positive_Matrix_With_Prefactor.hxx"

#include <boost/filesystem.hpp>
#include <algorithm>

void write_output_files(
  const boost::filesystem::path &outdir,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  auto max_normalization(
    std::max_element(normalization.begin(), normalization.end()));
  int64_t max_index(std::distance(normalization.begin(), max_normalization));
  El::BigFloat objective_const(objectives.at(max_index)
                               / normalization.at(max_index));
  std::vector<El::BigFloat> dual_objective_b;
  dual_objective_b.reserve(normalization.size() - 1);
  for(size_t index = 0; index < normalization.size(); ++index)
    {
      if(index != max_index)
        {
          dual_objective_b.push_back(objectives.at(index)
                                     - normalization.at(index) * objective_const);
        }
    }

  std::cout << "Objective: "
            << objective_const << "\n";
  for(auto &o : dual_objective_b)
    std::cout << " " << o << "\n";
}

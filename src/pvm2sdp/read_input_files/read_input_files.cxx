#include "../../sdp_convert.hxx"

#include <boost/filesystem.hpp>

std::vector<boost::filesystem::path>
read_file_list(const boost::filesystem::path &input_file);

void read_xml_input(const boost::filesystem::path &input_file,
                    El::BigFloat &objective_const,
                    std::vector<El::BigFloat> &dual_objectives_b,
                    std::vector<Dual_Constraint_Group> &dual_constraint_groups,
                    std::vector<size_t> &indices, size_t &num_processed);

void read_input_files(
  const std::vector<boost::filesystem::path> &input_files,
  El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  std::vector<size_t> &indices, size_t &num_processed);

void read_input_files(
  const std::vector<boost::filesystem::path> &input_files,
  El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  std::vector<size_t> &indices)
{
  size_t num_processed(0);

  read_input_files(input_files, objective_const, dual_objectives_b,
                   dual_constraint_groups, indices, num_processed);
}

void read_input_files(
  const std::vector<boost::filesystem::path> &input_files,
  El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  std::vector<size_t> &indices, size_t &num_processed)
{
  for(auto &input_file : input_files)
    {
      if(input_file.empty())
        {
          continue;
        }
      if(input_file.extension() == ".nsv")
        {
          read_input_files(read_file_list(input_file), objective_const,
                           dual_objectives_b, dual_constraint_groups, indices,
                           num_processed);
        }
      else
        {
          read_xml_input(input_file, objective_const, dual_objectives_b,
                         dual_constraint_groups, indices, num_processed);
        }
    }
}

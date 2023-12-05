#include "sdp_convert/sdp_convert.hxx"
#include "sdp_read/sdp_read.hxx"

namespace fs = std::filesystem;

void read_xml_input(const fs::path &input_file, El::BigFloat &objective_const,
                    std::vector<El::BigFloat> &dual_objectives_b,
                    std::vector<Dual_Constraint_Group> &dual_constraint_groups,
                    size_t &num_processed);

void read_input_files(
  const std::vector<fs::path> &input_files, El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  size_t &num_processed);

void read_input_files(
  const std::vector<fs::path> &input_files, El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  size_t &num_processed)
{
  for(auto &input_file : input_files)
    {
      if(input_file.empty())
        {
          continue;
        }
      if(input_file.extension() == ".nsv")
        {
          read_input_files(read_nsv_file_list(input_file), objective_const,
                           dual_objectives_b, dual_constraint_groups,
                           num_processed);
        }
      else if(input_file.extension() == ".xml")
        {
          read_xml_input(input_file, objective_const, dual_objectives_b,
                         dual_constraint_groups, num_processed);
        }
      else
        {
          El::RuntimeError("Cannot parse input file: ", input_file,
                           ". Expected .nsv or .xml extension.");
        }
    }
}

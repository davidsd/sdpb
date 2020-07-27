#include "../Functional.hxx"

#include <boost/filesystem/fstream.hpp>

Functional::Functional(const boost::filesystem::path &polynomials_path,
                       const boost::filesystem::path &poles_path)
    : has_prefactor(!poles_path.empty())
{
  {
    boost::filesystem::ifstream polynomials_file(polynomials_path);
    size_t block_index, poly_index, degree_index;
    std::string coefficient_string;
    polynomials_file >> block_index >> poly_index >> degree_index
      >> coefficient_string;
    while(polynomials_file)
      {
        if(block_index >= blocks.size())
          {
            blocks.resize(block_index + 1);
          }
        blocks[block_index].assign_poly(poly_index, degree_index,
                                        coefficient_string);
        polynomials_file >> block_index >> poly_index >> degree_index
          >> coefficient_string;
      }
  }

  if(has_prefactor)
    {
      boost::filesystem::ifstream poles_file(poles_path);
      size_t block_index, pole_index;
      std::string coefficient_string;
      poles_file >> block_index >> coefficient_string;
      while(poles_file)
        {
          if(block_index >= blocks.size())
            {
              blocks.resize(block_index + 1);
            }
          blocks[block_index].poles.emplace_back(coefficient_string);
          poles_file >> block_index >> coefficient_string;
        }
    }
}

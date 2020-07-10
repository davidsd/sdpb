#include "../Functional.hxx"

#include <boost/filesystem/fstream.hpp>

Functional::Functional(const boost::filesystem::path &polynomials_path,
                       const boost::filesystem::path &poles_path,
                       const bool &Has_prefactor)
    : has_prefactor(Has_prefactor)
{
  boost::filesystem::ifstream polynomials_file(polynomials_path);
  int64_t spin, index, degree;
  std::string coefficient_string;
  polynomials_file >> spin >> index >> degree >> coefficient_string;
  while(polynomials_file)
    {
      const int64_t spin_offset(spin / 2), index_offset(index - 1),
        degree_offset(degree - 1);
      if(spin_offset >= blocks.size())
        {
          blocks.resize(spin_offset + 1);
        }
      if(blocks[spin_offset].p.size() <= index_offset)
        {
          blocks[spin_offset].p.resize(index_offset + 1);
        }
      if(blocks[spin_offset].p[index_offset].coefficients.size()
         <= degree_offset)
        {
          blocks[spin_offset].p[index_offset].coefficients.resize(
            degree_offset + 1, El::BigFloat(0));
        }
      blocks[spin_offset].p[index_offset].coefficients[degree_offset]
        = El::BigFloat(coefficient_string);

      polynomials_file >> spin >> index >> degree >> coefficient_string;
    }

  // TODO poles
}

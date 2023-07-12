#pragma once

#include <El.hpp>

#include <boost/filesystem/fstream.hpp>

inline void write_distmatrix(const El::DistMatrix<El::BigFloat> &matrix,
                             const boost::filesystem::path &path)
{
  boost::filesystem::ofstream stream;
  if(matrix.DistRank() == matrix.Root())
    {
      boost::filesystem::create_directories(path.parent_path());
      stream.open(path);
    }
  El::Print(matrix,
            std::to_string(matrix.Height()) + " "
              + std::to_string(matrix.Width()),
            "\n", stream);
  if(matrix.DistRank() == matrix.Root())
    {
      stream << "\n";
      if(!stream.good())
        {
          throw std::runtime_error("Error when writing to: " + path.string());
        }
    }
}

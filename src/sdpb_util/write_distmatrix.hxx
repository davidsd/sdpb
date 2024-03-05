#pragma once

#include <El.hpp>

inline void write_distmatrix(const El::DistMatrix<El::BigFloat> &matrix,
                             const std::filesystem::path &path)
{
  std::ofstream stream;
  if(matrix.DistRank() == matrix.Root())
    {
      std::filesystem::create_directories(path.parent_path());
      stream.open(path);
    }
  El::Print(matrix,
            std::to_string(matrix.Height()) + " "
              + std::to_string(matrix.Width()),
            "\n", stream);
  if(matrix.DistRank() == matrix.Root())
    {
      stream << "\n";
      ASSERT(stream.good(), "Error when writing to: ", path);
    }
}

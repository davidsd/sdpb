#pragma once

#include "Block.hxx"

#include <boost/filesystem.hpp>

struct Functional
{
  std::vector<Block> blocks;
  bool has_prefactor;

  Functional(const boost::filesystem::path &polynomials_path);
  Functional(const boost::filesystem::path &polynomials_path,
             const boost::filesystem::path &poles_path);
  std::vector<El::BigFloat>
  eval(const std::vector<El::BigFloat> &coords,
       const std::vector<std::vector<El::BigFloat>> &optimals);

  El::BigFloat prefactor(const El::BigFloat &x);
};

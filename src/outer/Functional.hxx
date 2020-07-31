#pragma once

#include "Block.hxx"
#include "../sdp2input/Boost_Float.hxx"

#include <boost/filesystem.hpp>

struct Functional
{
  std::vector<Block> blocks;
  Boost_Float prefactor_power;

  Functional(const boost::filesystem::path &polynomials_path,
             const boost::filesystem::path &prefactor_path);
  Functional(const boost::filesystem::path &polynomials_path,
             const boost::filesystem::path &prefactor_path,
             const boost::filesystem::path &poles_path);
  std::vector<El::BigFloat>
  eval(const std::vector<El::BigFloat> &coords,
       const std::vector<std::vector<El::BigFloat>> &optimals) const;

  El::BigFloat prefactor(const El::BigFloat &x) const;
};

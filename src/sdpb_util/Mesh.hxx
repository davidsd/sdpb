#pragma once

#include <El.hpp>

#include <array>
#include <memory>
#include <iostream>

struct Mesh
{
  std::array<El::BigFloat, 5> x, f;
  std::unique_ptr<Mesh> lower, upper;

  Mesh(const El::BigFloat &x_0, const El::BigFloat &x_2,
       const El::BigFloat &x_4, const El::BigFloat &f_0,
       const El::BigFloat &f_2, const El::BigFloat &f_4,
       const std::function<El::BigFloat(const El::BigFloat &x)> &fn,
       const El::BigFloat &mesh_threshold,
       const El::BigFloat &block_epsilon);
  
  Mesh(const El::BigFloat &x_0, const El::BigFloat &x_4,
       const std::function<El::BigFloat(const El::BigFloat &x)> &fn,
       const El::BigFloat &mesh_threshold,
       const El::BigFloat &block_epsilon)
      : Mesh(x_0, (x_0 + x_4) / 2, x_4, fn(x_0), fn((x_0 + x_4) / 2), fn(x_4),
             fn, mesh_threshold, block_epsilon)
  {}
};

std::ostream & operator<<(std::ostream &os, const Mesh &mesh);

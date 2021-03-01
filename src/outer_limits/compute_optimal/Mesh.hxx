#pragma once

#include <El.hpp>

#include <array>
#include <memory>
#include <iostream>

struct Mesh
{
  std::array<El::BigFloat, 5> x, f;
  std::unique_ptr<Mesh> lower, upper;

  Mesh(const El::BigFloat &x_0, const El::BigFloat &x_3,
       const El::BigFloat &x_5, const El::BigFloat &f_0,
       const El::BigFloat &f_3, const El::BigFloat &f_5,
       const std::function<El::BigFloat(const El::BigFloat &x)> &fn,
       const El::BigFloat &epsilon,
       const El::BigFloat &block_epsilon);
  
  Mesh(const El::BigFloat &x_0, const El::BigFloat &x_5,
       const std::function<El::BigFloat(const El::BigFloat &x)> &fn,
       const El::BigFloat &epsilon,
       const El::BigFloat &block_epsilon)
      : Mesh(x_0, (x_0 + x_5) / 2, x_5, fn(x_0), fn((x_0 + x_5) / 2), fn(x_5),
             fn, epsilon, block_epsilon)
  {}
};

std::ostream & operator<<(std::ostream &os, const Mesh &mesh);

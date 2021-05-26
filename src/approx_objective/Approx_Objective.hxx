#pragma once

#include <El.hpp>

struct Approx_Objective
{
  El::BigFloat objective, dobjective, ddobjective;
  Approx_Objective(const El::BigFloat &Objective,
                   const El::BigFloat &Dobjective,
                   const El::BigFloat &Ddobjective)
      : objective(Objective), dobjective(Dobjective), ddobjective(Ddobjective)
  {}
};

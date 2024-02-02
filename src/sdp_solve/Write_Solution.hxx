#pragma once

#include <iostream>
#include <string>

struct Write_Solution
{
  std::string input_string;
  bool matrix_X=false, matrix_Y=false, vector_x=false, vector_y=false, vector_z=false;

  Write_Solution(const std::string &input);
  Write_Solution()=default;
};

inline std::ostream &
operator<<(std::ostream &os, const Write_Solution &write_solution)
{
  return os << write_solution.input_string;
}

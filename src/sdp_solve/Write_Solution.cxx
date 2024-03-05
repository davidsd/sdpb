#include "Write_Solution.hxx"

#include "sdpb_util/assert.hxx"

#include <boost/algorithm/string.hpp>

#include <stdexcept>
#include <vector>

Write_Solution::Write_Solution(const std::string &input) : input_string(input)
{
  std::vector<std::string> solutions;
  using namespace std::string_literals;
  boost::split(solutions, input_string, boost::is_any_of(", "s));
  for(auto &solution : solutions)
    {
      switch(solution.size())
        {
        case 0: break;
        case 1:
          switch(solution.front())
            {
            case 'x': vector_x = true; break;
            case 'y': vector_y = true; break;
            case 'z': vector_z = true; break;
            case 'X': matrix_X = true; break;
            case 'Y': matrix_Y = true; break;
            default:
              RUNTIME_ERROR(
                "Invalid argument for writeSolution.  Expected a comma "
                "separated list containing x, y, z, X, and/or Y, but found: "
                + solution);
            }
          break;
        default:
          RUNTIME_ERROR(
            "Invalid argument for writeSolution.  Expected a comma "
            "separated list containing x, y, X, and/or Y, but found: "
            + solution);
          break;
        }
    }
}

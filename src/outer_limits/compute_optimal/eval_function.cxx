#include <El.hpp>

#include <map>

El::BigFloat
eval_function(const El::BigFloat &infinity,
     const std::map<El::BigFloat, El::BigFloat> &f, const El::BigFloat &x)
{
  // lower_bound gives first element >=.  So in this case it is more
  // of the upper bound.
  auto next(f.lower_bound(x));
  if(next == f.end() || (next->first == infinity && x != infinity))
    {
      std::stringstream ss;
      ss << "Could not find this point in the function: '" << x << "'";
      throw std::runtime_error(ss.str());
    }
  if(next == f.begin())
    {
      // std::cout << "eval 0: "
      //           << x << " "
      //           << next->second
      //           << "\n";
      return next->second;
    }
  else
    {
      auto previous(next);
      --previous;
      const El::BigFloat delta(next->first - previous->first),
        dx((x - previous->first) / delta);

      // std::cout << "eval: "
      //           << x << " "
      //           << (1 - dx) * previous->second + dx * next->second
      //           << "\n";

      return (1 - dx) * previous->second + dx * next->second;
    }
}

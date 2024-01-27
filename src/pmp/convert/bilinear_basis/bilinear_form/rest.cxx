#include "Derivative_Term.hxx"
#include "../factorial.hxx"
#include "../accumulate_over_others.hxx"

#include <set>

std::set<Derivative_Term> dExp(const int64_t &k);

namespace
{
  Boost_Float
  log_rest(const int64_t &m, const Boost_Float &p,
           const std::vector<Boost_Float> &sorted_poles,
           const std::pair<std::vector<Boost_Float>::const_iterator,
                           std::vector<Boost_Float>::const_iterator> &equal_range,
           const int64_t &order)
  {
    return factorial(order - 1) * (order % 2 == 0 ? -1 : 1)
           * (m * pow(p, -order)
              - accumulate_over_others(
                  sorted_poles, equal_range, Boost_Float(0),
                  [&](const Boost_Float &sum, const Boost_Float &q) {
                    return sum + pow(p - q, -order);
                  }));
  }
}

Boost_Float
rest(const int64_t &m, const Boost_Float &p,
     const std::vector<Boost_Float> &sorted_poles,
     const std::pair<std::vector<Boost_Float>::const_iterator,
                     std::vector<Boost_Float>::const_iterator> &equal_range,
     const int64_t &k)
{
  Boost_Float result(0);
  for(auto &term : dExp(k))
    {
      Boost_Float product(term.constant);
      for(auto &power : term.powers)
        {
          product
            *= pow(log_rest(m, p, sorted_poles, equal_range, power.first),
                   power.second);
        }
      result += product;
    }
  return result / factorial(k);
}

#include "Derivative_Term.hxx"

#include <algorithm>
#include <set>

std::set<Derivative_Term> operator+(const std::set<Derivative_Term> &a,
                                    const std::set<Derivative_Term> &b)
{
  std::set<Derivative_Term> result;
  auto iter_a(a.begin()), iter_b(b.begin());
  while(iter_a != a.end())
    {
      if(iter_a->powers == iter_b->powers)
        {
          result.emplace(iter_a->constant + iter_b->constant, iter_a->powers);
          ++iter_a;
          ++iter_b;
        }
      else if(iter_a->powers < iter_b->powers)
        {
          result.insert(*iter_a);
          ++iter_a;
        }
      else
        {
          result.insert(*iter_b);
          ++iter_b;
        }
    }
  std::copy(iter_b, b.end(), std::inserter(result, result.begin()));
  return result;
}

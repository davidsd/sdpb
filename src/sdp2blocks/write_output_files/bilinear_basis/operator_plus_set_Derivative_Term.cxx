#include "Derivative_Term.hxx"

#include <algorithm>
#include <set>

void operator+=(std::set<Derivative_Term> &a,
                const std::set<Derivative_Term> &b)
{
  auto iter_a(a.begin()), iter_b(b.begin());
  while(iter_a != a.end())
    {
      if(iter_a->powers == iter_b->powers)
        {
          auto old_iter_a(iter_a);
          ++iter_a;
          a.emplace_hint(old_iter_a, old_iter_a->constant + iter_b->constant, iter_b->powers);
          a.erase(old_iter_a);
          ++iter_b;
        }
      else if(iter_a->powers < iter_b->powers)
        {
          ++iter_a;
        }
      else
        {
          a.insert(*iter_b);
          ++iter_b;
        }
    }
  std::copy(iter_b, b.end(), std::inserter(a, a.end()));
}

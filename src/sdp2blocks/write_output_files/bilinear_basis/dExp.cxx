#include "Derivative_Term.hxx"

std::set<Derivative_Term> dExp(const int64_t &k)
{
  std::set<Derivative_Term> result;
  result.emplace(1,std::map<int64_t,int64_t>());

  for(int64_t iteration=0; iteration<k; ++iteration)
    {
      std::set<Derivative_Term> new_terms;
      for(auto &term: result)
        {
          new_terms=new_terms+term.derivative();
          // std::swap(new_terms,new_terms+term.derivative());
        }
      std::swap(new_terms,result);
    }
  return result;
}

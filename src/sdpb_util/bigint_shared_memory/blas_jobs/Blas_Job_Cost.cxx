#include "Blas_Job_Cost.hxx"

#include <tuple>

Blas_Job_Cost::Blas_Job_Cost() : cost(0), blas_calls(0) {}

Blas_Job_Cost::Blas_Job_Cost(size_t elements, size_t blas_calls)
    : cost(elements), blas_calls(blas_calls)
{}

bool Blas_Job_Cost::operator<(const Blas_Job_Cost &rhs) const
{
  return std::tie(cost, blas_calls) < std::tie(rhs.cost, rhs.blas_calls);
}

Blas_Job_Cost Blas_Job_Cost::operator+(const Blas_Job_Cost &other) const
{
  return {cost + other.cost, blas_calls + other.blas_calls};
}

Blas_Job_Cost &Blas_Job_Cost::operator+=(const Blas_Job_Cost &other)
{
  cost += other.cost;
  blas_calls += other.blas_calls;
  return *this;
}

std::ostream &operator<<(std::ostream &os, const Blas_Job_Cost &cost)
{
  os << "cost: " << cost.cost << " blas_calls: " << cost.blas_calls;
  return os;
}

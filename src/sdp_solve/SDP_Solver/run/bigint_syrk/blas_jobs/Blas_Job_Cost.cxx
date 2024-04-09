#include "Blas_Job.hxx"

using Cost = Blas_Job::Cost;

Cost::Cost() : elements(0), blas_calls(0) {}

Cost::Cost(size_t elements, size_t blas_calls)
    : elements(elements), blas_calls(blas_calls)
{}

bool Cost::operator<(const Cost &rhs) const
{
  return std::tie(elements, blas_calls)
         < std::tie(rhs.elements, rhs.blas_calls);
}

Cost Cost::operator+(const Cost &other) const
{
  return {elements + other.elements, blas_calls + other.blas_calls};
}

Cost &Cost::operator+=(const Cost &other)
{
  elements += other.elements;
  blas_calls += other.blas_calls;
  return *this;
}

std::ostream &operator<<(std::ostream &os, const Cost &cost)
{
  os << "elements: " << cost.elements << " blas_calls: " << cost.blas_calls;
  return os;
}

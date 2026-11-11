#pragma once

#include <iostream>

struct Blas_Job_Cost
{
size_t cost;
size_t blas_calls;
Blas_Job_Cost();
Blas_Job_Cost(size_t elements, size_t blas_calls);
bool operator<(const Blas_Job_Cost &rhs) const;
Blas_Job_Cost operator+(Blas_Job_Cost const &other) const;
Blas_Job_Cost &operator+=(const Blas_Job_Cost &other);
friend std::ostream &operator<<(std::ostream &os, const Blas_Job_Cost &cost);
};

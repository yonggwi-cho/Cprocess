#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>

#define CHECK(cond)                                                       \
  do {                                                                    \
    if (!(cond)) {                                                        \
      std::fprintf(stderr, "FAILED %s:%d: %s\n", __FILE__, __LINE__, #cond); \
      std::exit(1);                                                       \
    }                                                                     \
  } while (0)

#define CHECK_NEAR(a, b, tol)                                              \
  do {                                                                     \
    const double a_ = (a), b_ = (b), t_ = (tol);                           \
    if (!(std::fabs(a_ - b_) <= t_)) {                                     \
      std::fprintf(stderr, "FAILED %s:%d: |%g - %g| = %g > %g\n", __FILE__, \
                   __LINE__, a_, b_, std::fabs(a_ - b_), t_);              \
      std::exit(1);                                                        \
    }                                                                      \
  } while (0)

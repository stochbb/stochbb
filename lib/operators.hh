#ifndef __SSB_OPERATORS_HH__
#define __SSB_OPERATORS_HH__

#include "randomvariable.hh"
#include "chain.hh"

inline sbb::Chain operator+(const sbb::RandomVariable &a, const sbb::RandomVariable &b) {
  return sbb::Chain(a,b);
}

#endif // __SSB_OPERATORS_HH__

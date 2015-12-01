#ifndef __SSB_OPERATORS_HH__
#define __SSB_OPERATORS_HH__

#include "api.hh"

/** Overloads the +-operator to construct a @c Chain.
 * @ingroup api */
inline sbb::Chain operator+(const sbb::Var &a, const sbb::Var &b) {
  return sbb::Chain(a,b);
}

namespace std {
/** Overloads the standard library @c max function to construct a @c Maximum RV.
 * @ingroup api */
inline  sbb::Minimum min(const sbb::Var &a, const sbb::Var &b) {
  return sbb::Minimum(a,b);
}

/** Overloads the standard library @c min function to construct a @c Minimum RV.
 * @ingroup api */
inline  sbb::Maximum max(const sbb::Var &a, const sbb::Var &b) {
  return sbb::Maximum(a,b);
}
}

#endif // __SSB_OPERATORS_HH__

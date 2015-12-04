#ifndef __SSB_OPERATORS_HH__
#define __SSB_OPERATORS_HH__

#include "api.hh"
#include <vector>

namespace sbb {

/** Constructs a random variable \f$Y\f$, as the minimum of the given random variables.
 * In contrast to the @c MinimumObj, it handles partial dependency between the random variables
 * if they are formed as sums of independent random variables by first separating the common
 * part. For example
 * \f[
 *  Y = min(X_1+X_2,X_1+X_3) \longrightarrow Y = X_1+min(X_2,X_3)\,,
 * \f]
 * where \f$X_1, X_2\f$ and \f$X_3\f$ are mutually independent random variables. */
Var minimum(const std::vector<Var> &variables);

/** Constructs a random variable \f$Y\f$, as the maximum of the given random variables.
 * In contrast to the @c MaximumObj, it handles partial dependency between the random variables
 * if they are formed as sums of independent random variables by first separating the common
 * part. For example
 * \f[
 *  Y = max(X_1+X_2,X_1+X_3) \longrightarrow Y = X_1+max(X_2,X_3)\,,
 * \f]
 * where \f$X_1, X_2\f$ and \f$X_3\f$ are mutually independent random variables. */
Var maximum(const std::vector<Var> &variables);

}

/** Overloads the +-operator to construct a @c Chain. */
inline sbb::Chain operator+(const sbb::Var &a, const sbb::Var &b) {
  std::vector<sbb::Var> args; args.reserve(2);
  args.push_back(a); args.push_back(b);
  return sbb::Chain(args);
}


namespace std {
/** Overloads the standard library @c max function to construct a @c Maximum RV. */
inline  sbb::Var min(const sbb::Var &a, const sbb::Var &b) {
  std::vector<sbb::Var> args; args.reserve(2); args.push_back(a); args.push_back(b);
  return sbb::minimum(args);
}

/** Overloads the standard library @c min function to construct a @c Minimum RV. */
inline  sbb::Var max(const sbb::Var &a, const sbb::Var &b) {
  std::vector<sbb::Var> args; args.reserve(2); args.push_back(a); args.push_back(b);
  return sbb::maximum(args);
}
}

#endif // __SSB_OPERATORS_HH__

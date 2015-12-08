#ifndef __SSB_OPERATORS_HH__
#define __SSB_OPERATORS_HH__

#include "api.hh"
#include <vector>

namespace sbb {

/** Tests if the vector of variables are mutually independent.
 * @param vars The vector of variables.
 * @returns @c true If the given variables are mutually independent. */
bool independent(const std::vector<Var> &vars);

/** Constructs a normal-distributed random variable for the given mean and standard deviation. */
Var normal(double mu, double sigma, const std::string &name="");
/** Constructs a compound-normal distributed random variable for the given mean and
 * standard deviation. Where the mean is a random variable too. */
Var normal(const Var &mu, double sigma, const std::string &name="");
/** Constructs a compound-normal distributed random variable for the given mean and
 * standard deviation. Where the standard deviation is a random variable too. */
Var normal(double mu, const Var &sigma, const std::string &name="");
/** Constructs a compound-normal distributed random variable for the given mean and
 * standard deviation. Where both, the mean and standard deviation are a random variables. */
Var normal(const Var& mu, const Var &sigma, const std::string &name="");

/** Constructs a gamma-distributed random variable for the given shape and scale parameters. */
Var gamma(double k, double theta, const std::string &name="");
/** Constructs a compound-gamma distributed random variable for the given shape and scale parameters.
 * Where the shape parameter is a random variable too. */
Var gamma(const Var &k, double theta, const std::string &name="");
/** Constructs a compound-gamma distributed random variable for the given shape and scale parameters.
 * Where the scale parameter is a random variable too. */
Var gamma(double k, const Var &theta, const std::string &name="");
/** Constructs a compound-gamma distributed random variable for the given shape and scale parameters.
 * Where both, the shape and scale parameter are random variables. */
Var gamma(const Var& k, const Var &theta, const std::string &name="");

/** Constructs a random variable \f$Y\f$ as the sum of the given random variables.
 * @param vars */
Var chain(const std::vector<Var> &vars);

/** Constructs a random variable \f$Y\f$ as the minimum of the given random variables.
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

/** Constructs an affine transform of the given variable. */
Var affine(const Var &var, double scale, double shift);
}

/** Overloads the +-operator to construct a @c Chain.
 * See also @c sbb::chain() function. */
inline sbb::Var operator +(const sbb::Var &a, const sbb::Var &b) {
  std::vector<sbb::Var> args; args.reserve(2);
  args.push_back(a); args.push_back(b);
  return sbb::chain(args);
}

/** Overloads the +-operator to construct an affine transform of the variable. */
inline sbb::Var operator +(const sbb::Var &var, double b) {
  return sbb::affine(var, 1, b);
}
/** Overloads the +-operator to construct an affine transform of the variable. */
inline sbb::Var operator +(double b, const sbb::Var &var) {
  return sbb::affine(var, 1, b);
}
/** Overloads the *-operator to construct an affine transform of the variable. */
inline sbb::Var operator *(double a, const sbb::Var &var) {
  return sbb::affine(var, a, 0);
}
/** Overloads the *-operator to construct an affine transform of the variable. */
inline sbb::Var operator *(const sbb::Var &var, double a) {
  return sbb::affine(var, a, 0);
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

#ifndef __SSB_OPERATORS_HH__
#define __SSB_OPERATORS_HH__

#include <ostream>
#include "api.hh"
#include <vector>

namespace stochbb {

/** Tests if the vector of variables are mutually independent.
 * @param vars The vector of variables.
 * @returns @c true If the given variables are mutually independent.
 * @ingroup api */
bool independent(const std::vector<Var> &vars);

/** Constructs a delta-distributed random variable (constant) located as the given value.
 * @ingroup api */
Var delta(double value);

/** Constructs a uniformly-distributed random variable on the interval [a,b].
 * @ingroup api */
Var uniform(double a, double b, const std::string &name="");

/** Constructs a normal-distributed random variable for the given mean and standard deviation.
 * @ingroup api */
Var normal(double mu, double sigma, const std::string &name="");
/** Constructs a compound-normal distributed random variable for the given mean and
 * standard deviation. Where the mean is a random variable too.
 * @ingroup api */
Var normal(const Var &mu, double sigma, const std::string &name="");
/** Constructs a compound-normal distributed random variable for the given mean and
 * standard deviation. Where the standard deviation is a random variable too.
 * @ingroup api */
Var normal(double mu, const Var &sigma, const std::string &name="");
/** Constructs a compound-normal distributed random variable for the given mean and
 * standard deviation. Where both, the mean and standard deviation are a random variables.
 * @ingroup api */
Var normal(const Var& mu, const Var &sigma, const std::string &name="");

/** Constructs a gamma-distributed random variable for the given shape and scale parameters.
 * @ingroup api */
Var gamma(double k, double theta, const std::string &name="");
/** Constructs a compound-gamma distributed random variable for the given shape and scale parameters.
 * Where the shape parameter is a random variable too.
 * @ingroup api */
Var gamma(const Var &k, double theta, const std::string &name="");
/** Constructs a compound-gamma distributed random variable for the given shape and scale parameters.
 * Where the scale parameter is a random variable too.
 * @ingroup api */
Var gamma(double k, const Var &theta, const std::string &name="");
/** Constructs a compound-gamma distributed random variable for the given shape and scale parameters.
 * Where both, the shape and scale parameter are random variables.
 * @ingroup api */
Var gamma(const Var& k, const Var &theta, const std::string &name="");

/** Constructs an inverse gamma-distributed random variable for the given shape and scale parameters.
 * @ingroup api */
Var invgamma(double alpha, double beta, const std::string &name="");
/** Constructs a compound inverse gamma distributed random variable for the given shape and scale
 * parameters. Where the shape parameter is a random variable too.
 * @ingroup api */
Var invgamma(const Var &alpha, double beta, const std::string &name="");
/** Constructs a compound inverse gamma distributed random variable for the given shape and scale
 * parameters. Where the scale parameter is a random variable too.
 * @ingroup api */
Var invgamma(double alpha, const Var &beta, const std::string &name="");
/** Constructs a compound inverse gamma distributed random variable for the given shape and scale
 * parameters. Where both, the shape and scale parameter are random variables.
 * @ingroup api */
Var invgamma(const Var& alpha, const Var &beta, const std::string &name="");

/** Constructs a Weibull-distributed random variable for the given shape and scale parameters.
 * @ingroup api */
Var weibull(double k, double lambda, const std::string &name="");
/** Constructs a compound Weibull-distributed random variable for the given shape and scale parameters.
 * @ingroup api */
Var weibull(const Var& k, double lambda, const std::string &name="");
/** Constructs a compound Weibull-distributed random variable for the given shape and scale parameters.
 * @ingroup api */
Var weibull(double k, const Var& lambda, const std::string &name="");
/** Constructs a compound Weibull-distributed random variable for the given shape and scale parameters.
 * @ingroup api */
Var weibull(const Var& k, const Var& lambda, const std::string &name="");

/** Constructs a random variable \f$Y\f$ as the sum of the given random variables.
 * That is
 * \f[
 *  Y = \sum_i X_i\,.
 * \f]
 *
 * @param vars Specifies the RVs to sum.
 * @ingroup api */
Var chain(const std::vector<Var> &vars);

/** Constructs a random variable \f$Y\f$ as the sum of the given random variables.
 * That is
 * \f[
 *  Y = X_1+X_2\,.
 * \f]
 * @ingroup api */
Var chain(const Var &X1, const Var &X2);

/** Constructs a random variable \f$Y\f$ as the sum of the given random variables.
 * That is
 * \f[
 *  Y = X_1+X_2+X_3\,.
 * \f]
 * @ingroup api */
Var chain(const Var &X1, const Var &X2, const Var &X3);

/** Constructs a random variable \f$Y\f$ as the minimum of the given random variables.
 * That is
 * \f[
 *  Y = min(X_1, X_2, ..., X_N)\,.
 * \f]
 *
 * In contrast to the direct instantiation of the @c Minimum class, it handles partial dependency
 * between the random variables if they are formed as sums of independent random variables by first
 * separating the common part. For example
 * \f[
 *  Y = min(X_1+X_2,X_1+X_3) \longrightarrow Y = X_1+min(X_2,X_3)\,,
 * \f]
 * where \f$X_1, X_2\f$ and \f$X_3\f$ are mutually independent random variables.
 * @ingroup api */
Var minimum(const std::vector<Var> &variables);

/** Constructs a random variable \f$Y\f$ as the minimum of the given random variables.
 * That is
 * \f[
 *  Y = min(X_1, X_2)\,.
 * \f]
 *
 * In contrast to the direct instantiation of the @c Minimum class, it handles partial dependency
 * between the random variables if they are formed as sums of independent random variables by first
 * separating the common part. For example
 * \f[
 *  Y = min(X_1+X_2,X_1+X_3) \longrightarrow Y = X_1+min(X_2,X_3)\,,
 * \f]
 * where \f$X_1, X_2\f$ and \f$X_3\f$ are mutually independent random variables.
 * @ingroup api */
Var minimum(const Var &X1, const Var &X2);

/** Constructs a random variable \f$Y\f$ as the minimum of the given random variables.
 * That is
 * \f[
 *  Y = min(X_1, X_2, X_3)\,.
 * \f]
 *
 * In contrast to the direct instantiation of the @c Minimum class, it handles partial dependency
 * between the random variables if they are formed as sums of independent random variables by first
 * separating the common part. For example
 * \f[
 *  Y = min(X_1+X_2,X_1+X_3) \longrightarrow Y = X_1+min(X_2,X_3)\,,
 * \f]
 * where \f$X_1, X_2\f$ and \f$X_3\f$ are mutually independent random variables.
 * @ingroup api */
Var minimum(const Var &X1, const Var &X2, const Var &X3);

/** Constructs a random variable \f$Y\f$, as the maximum of the given random variables.
 * That is
 * \f[
 *  Y = max(X_1, X_2, ..., X_N)\,.
 * \f]
 *
 * In contrast to the direct instantiation of the @c Maximum class, it handles partial dependency
 * between the random variables if they are formed as sums of independent random variables by first
 * separating the common part. For example
 * \f[
 *  Y = max(X_1+X_2,X_1+X_3) \longrightarrow Y = X_1+max(X_2,X_3)\,,
 * \f]
 * where \f$X_1, X_2\f$ and \f$X_3\f$ are mutually independent random variables.
 * @ingroup api */
Var maximum(const std::vector<Var> &variables);

/** Constructs a random variable \f$Y\f$, as the maximum of the given random variables.
 * That is
 * \f[
 *  Y = max(X_1, X_2)\,.
 * \f]
 *
 * In contrast to the direct instantiation of the @c Maximum class, it handles partial dependency
 * between the random variables if they are formed as sums of independent random variables by first
 * separating the common part. For example
 * \f[
 *  Y = max(X_1+X_2,X_1+X_3) \longrightarrow Y = X_1+max(X_2,X_3)\,,
 * \f]
 * where \f$X_1, X_2\f$ and \f$X_3\f$ are mutually independent random variables.
 * @ingroup api */
Var maximum(const Var &X1, const Var &X2);

/** Constructs a random variable \f$Y\f$, as the maximum of the given random variables.
 * That is
 * \f[
 *  Y = max(X_1, X_2, X_3)\,.
 * \f]
 *
 * In contrast to the direct instantiation of the @c Maximum class, it handles partial dependency
 * between the random variables if they are formed as sums of independent random variables by first
 * separating the common part. For example
 * \f[
 *  Y = max(X_1+X_2,X_1+X_3) \longrightarrow Y = X_1+max(X_2,X_3)\,,
 * \f]
 * where \f$X_1, X_2\f$ and \f$X_3\f$ are mutually independent random variables.
 * @ingroup api */
Var maximum(const Var &X1, const Var &X2, const Var &X3);

/** Constructs an affine transformed of the given variable.
 * @ingroup api */
Var affine(const Var &var, double scale, double shift);

/** Constructs a mixture of the given variables with assiciated weights.
 * @ingroup api */
Var mixture(double wX1, const Var &X1, double wX2, const Var &X2);

/** Constructs a mixture of the given variables with assiciated weights.
 * @ingroup api */
Var mixture(double wX1, const Var &X1, double wX2, const Var &X2, double wX3, const Var &X3);

/** Constructs a mixture of the given variables with assiciated weights.
 * @ingroup api */
Var mixture(const std::vector<double> &weights, const std::vector<Var> &variables);

/** Constructs a conditional random variable. That is
 * \f[
 *  Z = \begin{cases}Y_1 & if\, X_1<X_2 \\ Y_2 & else.\end{cases}
 * \f]
 * @ingroup api */
Var conditional(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2);

/** Constructs a conditional chained random variable. That is
 * \f[
 *  Z = \begin{cases}X_1+Y_1 & if\, X_1<X_2 \\ X_2+Y_2 & else.\end{cases}
 * \f]
 * @ingroup api */
Var condchain(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2);
}

/** Overloads the +-operator to construct a @c Chain.
 * See also @c sbb::chain() function. */
inline stochbb::Var operator +(const stochbb::Var &a, const stochbb::Var &b) {
  std::vector<stochbb::Var> args; args.reserve(2);
  args.push_back(a); args.push_back(b);
  return stochbb::chain(args);
}

/** Overloads the +-operator to construct an affine transform of the variable. */
inline stochbb::Var operator +(const stochbb::Var &var, double b) {
  return stochbb::affine(var, 1, b);
}
/** Overloads the +-operator to construct an affine transform of the variable. */
inline stochbb::Var operator +(double b, const stochbb::Var &var) {
  return stochbb::affine(var, 1, b);
}
/** Overloads the *-operator to construct an affine transform of the variable. */
inline stochbb::Var operator *(double a, const stochbb::Var &var) {
  return stochbb::affine(var, a, 0);
}
/** Overloads the *-operator to construct an affine transform of the variable. */
inline stochbb::Var operator *(const stochbb::Var &var, double a) {
  return stochbb::affine(var, a, 0);
}

inline std::ostream &operator<<(std::ostream &stream, const stochbb::Container &x) {
  if (x.isNull()) {
    stream << "<null>";
  } else {
    x->print(stream);
  }
  return stream;
}

namespace std {

/** Overloads the standard library @c max function to construct a @c Maximum RV. */
inline  stochbb::Var min(const stochbb::Var &a, const stochbb::Var &b) {
  std::vector<stochbb::Var> args; args.reserve(2); args.push_back(a); args.push_back(b);
  return stochbb::minimum(args);
}

/** Overloads the standard library @c min function to construct a @c Minimum RV. */
inline  stochbb::Var max(const stochbb::Var &a, const stochbb::Var &b) {
  std::vector<stochbb::Var> args; args.reserve(2); args.push_back(a); args.push_back(b);
  return stochbb::maximum(args);
}

}


#endif // __SSB_OPERATORS_HH__

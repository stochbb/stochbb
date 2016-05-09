/** @defgroup math Core mathematical functions
 * @ingroup internal */

#ifndef __SBB_MATH_HH__
#define __SBB_MATH_HH__

#include <cmath>

namespace stochbb {

/** The regularized lower-incomplete gamma function.
 * @ingroup math */
double gamma_li(double s, double z);
/** The regularized upper-incomplete gamma function.
 * @ingroup math */
double gamma_ui(double s, double z);
/** Logarithm of the gamma function.
 * @ingroup math */
double lgamma(double z);

/** The error function.
 * @ingroup math */
double erfc(double x);
/** The inverse error function.
 * @ingroup math */
double erfinv(double p);

/** Regularized incomplete beta function.
 * @ingroup math */
double betai(double a, double b, double x);

/** Inverse regularized incomplete beta function.
 * @ingroup math */
double invbetai(double a, double b, double x);

/** PDF of the standard normal distribution.
 * @ingroup math */
double dnorm(double x);
/** CDF of the standard normal distribution.
 * @ingroup math */
double pnorm(double x);
/** Inverse CDF of the standard normal distribution (quantile function).
 * @ingroup math */
double qnorm(double p);

/** PDF of the gamma distribution.
 * @ingroup math */
double dgamma(double x, double k, double theta);
/** CDF of the gamma distribution.
 * @ingroup math */
double pgamma(double x, double k, double theta);
/** Inverse CDF of the gamma distribution (quantile function).
 * @ingroup math */
double qgamma(double p, double k, double theta);

/** PDF of the inverse gamma distribution.
 * @ingroup math */
double dinvgamma(double x, double alpha, double beta);
/** CDF of the inverse gamma distribution.
 * @ingroup math */
double pinvgamma(double x, double alpha, double beta);
/** Inverse CDF of the inverse gamma distribution (quantile function).
 * @ingroup math */
double qinvgamma(double p, double alpha, double beta);

/** PDF of the Weibull distribution.
 * @ingroup math */
double dweibull(double x, double k, double lambda);
/** CDF of the Weibull distribution.
 * @ingroup math */
double pweibull(double x, double k, double lambda);
/** Inverse CDF of the Weibull distribution (quantile function).
 * @ingroup math */
double qweibull(double p, double k, double lambda);

}

#endif // __SBB_MATH_HH__

#ifndef __SBB_MATH_HH__
#define __SBB_MATH_HH__

#include <cmath>

namespace stochbb {

/** The regularized lower-incomplete gamma function. */
double gamma_li(double s, double z);
/** The regularized upper-incomplete gamma function. */
double gamma_ui(double s, double z);
/** Logarithm of the gamma function. */
double lgamma(double z);

/** The error function. */
double erfc(double x);
/** The inverse error function. */
double erfinv(double p);

/** Regularized incomplete beta function. */
double betai(double a, double b, double x);

/** PDF of the standard normal distribution. */
double dnorm(double x);
/** CDF of the standard normal distribution. */
double pnorm(double x);
/** Inverse CDF of the standard normal distribution (quantile function). */
double qnorm(double p);

/** PDF of the gamma distribution. */
double dgamma(double x, double k, double theta);
/** CDF of the gamma distribution. */
double pgamma(double x, double k, double theta);
/** Inverse CDF of the gamma distribution (quantile function). */
double qgamma(double p, double k, double theta);

/** PDF of the inverse gamma distribution. */
double dinvgamma(double x, double alpha, double beta);
/** CDF of the inverse gamma distribution. */
double pinvgamma(double x, double alpha, double beta);
/** Inverse CDF of the inverse gamma distribution (quantile function). */
double qinvgamma(double p, double alpha, double beta);

/** PDF of the Weibull distribution. */
double dweibull(double x, double k, double lambda);
/** CDF of the Weibull distribution. */
double pweibull(double x, double k, double lambda);
/** Inverse CDF of the Weibull distribution (quantile function). */
double qweibull(double p, double k, double lambda);

}

#endif // __SBB_MATH_HH__

#include "math.hh"
#include <limits>
#include <cstddef>

using namespace stochbb;

/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 */
double
stochbb::lgamma(double z) {
  double x = 0;
  x += 0.1659470187408462e-06 / (z+7);
  x += 0.9934937113930748e-05 / (z+6);
  x -= 0.1385710331296526     / (z+5);
  x += 12.50734324009056      / (z+4);
  x -= 176.6150291498386      / (z+3);
  x += 771.3234287757674      / (z+2);
  x -= 1259.139216722289      / (z+1);
  x += 676.5203681218835      / z;
  x += 0.9999999999995183;
  return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

/* complementary error function
 * \frac{2}{\sqrt{\pi}} \int_x^{\infty} e^{-t^2} dt
 * AS66, 2nd algorithm, http://lib.stat.cmu.edu/apstat/66
 */
double
stochbb::erfc(double x) {
  const double p0 = 220.2068679123761;
  const double p1 = 221.2135961699311;
  const double p2 = 112.0792914978709;
  const double p3 = 33.912866078383;
  const double p4 = 6.37396220353165;
  const double p5 = .7003830644436881;
  const double p6 = .03526249659989109;
  const double q0 = 440.4137358247522;
  const double q1 = 793.8265125199484;
  const double q2 = 637.3336333788311;
  const double q3 = 296.5642487796737;
  const double q4 = 86.78073220294608;
  const double q5 = 16.06417757920695;
  const double q6 = 1.755667163182642;
  const double q7 = .08838834764831844;
  double expntl, z, p;
  z = fabs(x) * M_SQRT2;
  if (z > 37.) return x > 0.? 0. : 2.;
  expntl = exp(z * z * - .5);
  if (z < 10. / M_SQRT2) // for small z
      p = expntl * ((((((p6 * z + p5) * z + p4) * z + p3) * z + p2) * z + p1) * z + p0)
      / (((((((q7 * z + q6) * z + q5) * z + q4) * z + q3) * z + q2) * z + q1) * z + q0);
  else p = expntl / 2.506628274631001 / (z + 1. / (z + 2. / (z + 3. / (z + 4. / (z + .65)))));
  return x > 0.? 2. * p : 2. * (1. - p);
}

/* The following computes regularized incomplete gamma functions.
 * Formulas are taken from Wiki, with additional input from Numerical
 * Recipes in C (for modified Lentz's algorithm) and AS245
 * (http://lib.stat.cmu.edu/apstat/245).
 *
 * A good online calculator is available at:
 *
 *   http://www.danielsoper.com/statcalc/calc23.aspx
 *
 * It calculates upper incomplete gamma function, which equals
 * kf_gammaq(s,z)*tgamma(s).
 */

#define SBB_GAMMA_EPS 1e-14
#define SBB_TINY 1e-290

// regularized lower incomplete gamma function, by series expansion
double
_loggammap(double s, double z) {
  double sum, x;
  int k;
  for (k = 1, sum = x = 1.; k < 100; ++k) {
    sum += (x *= z / (s + k));
    if (x / sum < SBB_GAMMA_EPS) break;
  }
  return s * log(z) - z - stochbb::lgamma(s + 1.) + log(sum);
}

inline double
_gammap(double s, double z) {
  return std::exp(_loggammap(s,z));
}

// regularized upper incomplete gamma function, by continued fraction
double
_loggammaq(double s, double z) {
  int j;
  double C, D, f;
  f = 1. + z - s; C = f; D = 0.;
  // Modified Lentz's algorithm for computing continued fraction
  // See Numerical Recipes in C, 2nd edition, section 5.2
  for (j = 1; j < 100; ++j) {
    double a = j * (s - j), b = (j<<1) + 1 + z - s, d;
    D = b + a * D;
    if (D < SBB_TINY) D = SBB_TINY;
    C = b + a / C;
    if (C < SBB_TINY) C = SBB_TINY;
    D = 1. / D;
    d = C * D;
    f *= d;
    if (fabs(d - 1.) < SBB_GAMMA_EPS) break;
  }
  return s * log(z) - z - stochbb::lgamma(s) - log(f);
}

inline double
_gammaq(double s, double z) {
  return std::exp(_loggammaq(s,z));
}

double
stochbb::gamma_li(double s, double z) {
  return z <= 1. || z < s ? _gammap(s, z) : 1. - _gammaq(s, z);
}

double
stochbb::gamma_ui(double s, double z) {
  return z <= 1. || z < s? 1. - _gammap(s, z) : _gammaq(s, z);
}

/* Regularized incomplete beta function. The method is taken from
 * Numerical Recipe in C, 2nd edition, section 6.4. The following web
 * page calculates the incomplete beta function, which equals
 * kf_betai(a,b,x) * gamma(a) * gamma(b) / gamma(a+b):
 *
 *   http://www.danielsoper.com/statcalc/calc36.aspx
 */
double
_betai_aux(double a, double b, double x) {
  double C, D, f;
  int j;
  if (x == 0.) return 0.;
  if (x == 1.) return 1.;
  f = 1.; C = f; D = 0.;
  // Modified Lentz's algorithm for computing continued fraction
  for (j = 1; j < 200; ++j) {
    double aa, d;
    int m = j>>1;
    aa = (j&1)? -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1))
      : m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m));
    D = 1. + aa * D;
    if (D < SBB_TINY) D = SBB_TINY;
    C = 1. + aa / C;
    if (C < SBB_TINY) C = SBB_TINY;
    D = 1. / D;
    d = C * D;
    f *= d;
    if (fabs(d - 1.) < SBB_GAMMA_EPS) break;
  }
  return exp(stochbb::lgamma(a+b) - stochbb::lgamma(a) - stochbb::lgamma(b) + a * log(x) + b * log(1.-x)) / a / f;
}

double
stochbb::betai(double a, double b, double x) {
  return x < (a + 1.) / (a + b + 2.) ? _betai_aux(a, b, x) : 1. - _betai_aux(b, a, 1. - x);
}

double
stochbb::invbetai(double a, double b, double p) {
  // Perform binary search for x
  double lower=0, upper=1, x=0;
  while ((upper-lower)>1e-9) {
    x = (upper-lower)/2;
    if (betai(a,b,x) > p) {
      upper=x;
    } else {
      lower=x;
    }
  }
  return 0;
}

double
stochbb::erfinv(double x) {
  // beware that the logarithm argument must be
  // commputed as (1.0 - x) * (1.0 + x),
  // it must NOT be simplified as 1.0 - x * x as this
  // would induce rounding errors near the boundaries +/-1
  double w = - std::log((1.0 - x) * (1.0 + x));
  double p;

  if (w < 6.25) {
    w -= 3.125;
    p =  -3.6444120640178196996e-21;
    p =   -1.685059138182016589e-19 + p * w;
    p =   1.2858480715256400167e-18 + p * w;
    p =    1.115787767802518096e-17 + p * w;
    p =   -1.333171662854620906e-16 + p * w;
    p =   2.0972767875968561637e-17 + p * w;
    p =   6.6376381343583238325e-15 + p * w;
    p =  -4.0545662729752068639e-14 + p * w;
    p =  -8.1519341976054721522e-14 + p * w;
    p =   2.6335093153082322977e-12 + p * w;
    p =  -1.2975133253453532498e-11 + p * w;
    p =  -5.4154120542946279317e-11 + p * w;
    p =    1.051212273321532285e-09 + p * w;
    p =  -4.1126339803469836976e-09 + p * w;
    p =  -2.9070369957882005086e-08 + p * w;
    p =   4.2347877827932403518e-07 + p * w;
    p =  -1.3654692000834678645e-06 + p * w;
    p =  -1.3882523362786468719e-05 + p * w;
    p =    0.0001867342080340571352 + p * w;
    p =  -0.00074070253416626697512 + p * w;
    p =   -0.0060336708714301490533 + p * w;
    p =      0.24015818242558961693 + p * w;
    p =       1.6536545626831027356 + p * w;
  } else if (w < 16.0) {
    w = std::sqrt(w) - 3.25;
    p =   2.2137376921775787049e-09;
    p =   9.0756561938885390979e-08 + p * w;
    p =  -2.7517406297064545428e-07 + p * w;
    p =   1.8239629214389227755e-08 + p * w;
    p =   1.5027403968909827627e-06 + p * w;
    p =   -4.013867526981545969e-06 + p * w;
    p =   2.9234449089955446044e-06 + p * w;
    p =   1.2475304481671778723e-05 + p * w;
    p =  -4.7318229009055733981e-05 + p * w;
    p =   6.8284851459573175448e-05 + p * w;
    p =   2.4031110387097893999e-05 + p * w;
    p =   -0.0003550375203628474796 + p * w;
    p =   0.00095328937973738049703 + p * w;
    p =   -0.0016882755560235047313 + p * w;
    p =    0.0024914420961078508066 + p * w;
    p =   -0.0037512085075692412107 + p * w;
    p =     0.005370914553590063617 + p * w;
    p =       1.0052589676941592334 + p * w;
    p =       3.0838856104922207635 + p * w;
  } else if (std::numeric_limits<double>::infinity() != std::abs(w)) {
    w = std::sqrt(w) - 5.0;
    p =  -2.7109920616438573243e-11;
    p =  -2.5556418169965252055e-10 + p * w;
    p =   1.5076572693500548083e-09 + p * w;
    p =  -3.7894654401267369937e-09 + p * w;
    p =   7.6157012080783393804e-09 + p * w;
    p =  -1.4960026627149240478e-08 + p * w;
    p =   2.9147953450901080826e-08 + p * w;
    p =  -6.7711997758452339498e-08 + p * w;
    p =   2.2900482228026654717e-07 + p * w;
    p =  -9.9298272942317002539e-07 + p * w;
    p =   4.5260625972231537039e-06 + p * w;
    p =  -1.9681778105531670567e-05 + p * w;
    p =   7.5995277030017761139e-05 + p * w;
    p =  -0.00021503011930044477347 + p * w;
    p =  -0.00013871931833623122026 + p * w;
    p =       1.0103004648645343977 + p * w;
    p =       4.8499064014085844221 + p * w;
  } else {
    // this branch does not appears in the original code, it
    // was added because the previous branch does not handle
    // x = +/-1 correctly. In this case, w is positive infinity
    // and as the first coefficient (-2.71e-11) is negative.
    // Once the first multiplication is done, p becomes negative
    // infinity and remains so throughout the polynomial evaluation.
    // So the branch above incorrectly returns negative infinity
    // instead of the correct positive infinity.
    p = std::numeric_limits<double>::infinity();
  }
  // done.
  return p * x;
}


double
stochbb::dnorm(double x) {
  return std::exp(-x*x/2)/std::sqrt(2*M_PI);
}

double
stochbb::pnorm(double x) {
  return 0.5*(1+std::erf(x/std::sqrt(2)));
}

double
stochbb::qnorm(double p) {
  return std::sqrt(2) * stochbb::erfinv(2*p-1);
}

double
stochbb::dgamma(double x, double k, double theta) {
  if (x<=0) {
    return 0;
  }
  if ((x==0) && (k==1)) {
    return std::exp(-std::lgamma(k) -k*std::log(theta));
  } else if ((x==0) && (k>1)) {
    return 0;
  }
  return std::exp((k-1)*std::log(x) - x/theta -std::lgamma(k) -k*std::log(theta));
}

double
stochbb::pgamma(double x, double k, double theta) {
  if (x<=0) { return 0; }
  return stochbb::gamma_li(k, x/theta);
}

double
stochbb::qgamma(double p, double k, double theta) {
  // Solve p = pgamma(x, k, theta) by simple Newton
  // initial guess
  double x = k*theta;
  // p(x)
  double p2 = stochbb::pgamma(x, k, theta);
  // do Newton
  while (std::abs(p-p2)>1e-8) {
    // update x
    x += (p-p2)/stochbb::dgamma(x, k, theta);
    // get new p(x)
    p2 = stochbb::pgamma(x, k, theta);
  }
  // done.
  return x;
}


double
stochbb::dinvgamma(double x, double alpha, double beta) {
  if (x<=0) { return 0; }
  return std::exp( alpha*std::log(beta) - std::lgamma(alpha) - (alpha+1)*std::log(x) - beta/x );
}

double
stochbb::pinvgamma(double x, double alpha, double beta) {
  if (x<=0) { return 0; }
  return stochbb::gamma_ui(alpha, beta/x);
}

double
stochbb::qinvgamma(double p, double alpha, double beta) {
  // Solve p = pgamma(x, k, theta) by simple Newton
  // Initial quess
  double x = beta/(alpha+1);
  // get p(x)
  double p2 = stochbb::pinvgamma(x, alpha, beta);
  while (std::abs(p-p2)>1e-8) {
    // update x
    x += (p-p2)/stochbb::dinvgamma(x, alpha, beta);
    // update p(x)
    p2 = stochbb::pinvgamma(x, alpha, beta);
  }
  // done.
  return x;
}


double
stochbb::dweibull(double x, double k, double lambda) {
  if (x<0) { return 0; }
  return std::exp(std::log(k)-std::log(lambda)
                  + (k-1)*(std::log(x)-std::log(lambda))
                  - std::pow(x/lambda, k));
}

double
stochbb::pweibull(double x, double k, double lambda) {
  if (x<0) { return 0; }
  return 1-std::exp(std::pow(x/lambda,k));
}

double
stochbb::qweibull(double p, double k, double lambda) {
  return lambda*std::pow(std::log(1-p), 1./k);
}

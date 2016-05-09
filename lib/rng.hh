#ifndef __SBB_RNG_HH__
#define __SBB_RNG_HH__

#include "mersennetwister.hh"

namespace stochbb {

/** Singleton class of the random number generator using the @c MersenneTwister.
 * The singleton instance is not directly accessible, but random numbers can be drawn using the
 * static methods.
 * @ingroup math */
class RNG: protected MersenneTwister
{
protected:
  /** Hidden constructor. */
  RNG();

public:
  /** Returns a uniformly distributed random value on the interval [0,1]. */
  static double unif();
  /** Returns a uniformly distributed random value on the interval [a,b]. */
  static double unif(double a, double b);
  /** Returns a stdandard normal distributed random value. */
  static double norm();
  /** Returns a normal distributed random value with mean mu and standard deviation sigma. */
  static double norm(double mu, double sigma);
  /** Returns a Gamma(k, theta) distributed random variable. */
  static double gamma(double k, double theta);
  /** Returns a InvGamma(alpha, beta) distributed random variable. */
  static double invgamma(double alpha, double beta);

protected:
  /** Factory function. */
  static RNG *get();
  /** Singleton instance. */
  static RNG *_instance;
};

}

#endif // __SBB_RNG_HH__

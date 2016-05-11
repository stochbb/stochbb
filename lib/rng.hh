#ifndef __SBB_RNG_HH__
#define __SBB_RNG_HH__

#include <random>

namespace stochbb {

/** Singleton class of the random number generator using the Mersenne twister (@c std::mt19937_64)
 * RNG. This class can be used in conjecture with the distributions defined in C++ std. libs
 * 'random'.
 * @ingroup math */
class RNG: public std::mt19937_64
{
protected:
  /** Hidden constructor. */
  RNG();

public:
  /** Factory function. */
  static RNG &get();

protected:
  /** Singleton instance. */
  static RNG *_instance;
};

}

namespace std {

/** Implements the inverse gamma distribution. */
template <typename RealType=double>
class inverse_gamma_distribution: public gamma_distribution<RealType>
{
public:
  /** Constructs a new inverse gamma sampler with the given parameters. */
  explicit inverse_gamma_distribution(RealType alpha = 1, RealType beta = 1)
    : gamma_distribution<RealType>(alpha, beta) {}

  /** Constructs a new inverse gamma sampler with the given parameters. */
  explicit inverse_gamma_distribution(const typename gamma_distribution<RealType>::param_type& p)
    : gamma_distribution<RealType>(p) {}

  /** Generates inverse gamma distributed random values using the given RNG. */
  template<class URNG> RealType operator()(URNG& g) {
    return 1./gamma_distribution<RealType>::operator() (g);
  }

  /** Generates inverse gamma distributed random values using the given RNG and parameters. */
  template<class URNG> RealType operator()(URNG& g, const typename gamma_distribution<RealType>::param_type& param) {
    return 1./gamma_distribution<RealType>::operator() (g, param);
  }
};

}

#endif // __SBB_RNG_HH__

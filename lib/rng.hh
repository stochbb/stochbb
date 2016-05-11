#ifndef __SBB_RNG_HH__
#define __SBB_RNG_HH__

#include <random>

namespace stochbb {

/** Singleton class of the random number generator using the Mersenne twister (@c std::mt19937_64)
 * RNG. This class can be used in conjecture with the distributions defined in C++ std. lib
 * <random>.
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
  explicit inverse_gamma_distribution(RealType alpha = 1, RealType beta = 1)
    : gamma_distribution<RealType>(alpha, beta) {}

  explicit inverse_gamma_distribution(const typename gamma_distribution<RealType>::param_type& p)
    : gamma_distribution<RealType>(p) {}

  /** Generates inverse gamma distributed random values using the given RNG. */
  template<class URNG> RealType operator()(URNG& g);
  /** Generates inverse gamma distributed random values using the given RNG and parameters. */
  template<class URNG> RealType operator()(URNG& g, const typename gamma_distribution<RealType>::param_type& parm);
};


/* ********************************************************************************************* *
 * Template implementation
 * ********************************************************************************************* */
template<class _RealType>
template<class _URNG>
_RealType
inverse_gamma_distribution<_RealType>::operator() (_URNG &__g, const typename gamma_distribution<_RealType>::param_type &__p)
{
  return 1./gamma_distribution<_RealType>::operator() (__g, __p);
}

template<class _RealType>
template<class _URNG>
_RealType
inverse_gamma_distribution<_RealType>::operator() (_URNG &__g)
{
  return 1./gamma_distribution<_RealType>::operator() (__g);
}

}
#endif // __SBB_RNG_HH__

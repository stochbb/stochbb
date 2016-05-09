#ifndef __SBB_DIRECTMARGINALSAMPLER_HH__
#define __SBB_DIRECTMARGINALSAMPLER_HH__

#include <Eigen/Eigen>
#include "api.hh"

namespace stochbb {

/** Implements a simple and possibly fast marginal sampler for a single variable.
 * This marginal sampler uses the numerical approximation of the derived CDF and
 * its inverse to implement a generic sampler for a specified variable.
 *
 * This sampler uses the so-called "golden rule" to construct an approximative sampler
 * from an piece-wise linear CDF: Assume that probability function of a random variable
 * \f$X\f$ is \f$F_X(x)\f$, then \f$Y = F_X(X)\f$ is uniformly distributed on the interval
 * \f$[0,1]\f$. Consquently will \f$F_X^{-1}(Y)\f$ be distributed like \f$X\f$.
 * @ingroup internal */
class MarginalSamplerObj : public Object
{
public:
  /** Constructs a direct marginal sampler.
   * @param variable Specifies the random variable to sample from.
   * @param Tmin Specifies the lower-bound of the interval on which the CDF is evaluated.
   * @param Tmax Specifies the upper-bound of the interval on which the CDF is evaluated.
   * @param steps Specifies the number of bins on which the CDF is evaluated. */
  MarginalSamplerObj(const Var &variable, double Tmin, double Tmax, size_t steps);
  virtual void mark();

  /** Samples from the variable. */
  void sample(Eigen::Ref<Eigen::VectorXd> out);

protected:
  /** The marginal to sample from. */
  VarObj *_variable;
  /** The lower-bound of the interval on which the CDF was evaluated. */
  double _Tmin;
  /** The upper-bound of the interval on which the CDF was evaluated. */
  double _Tmax;
  /** The evaluated CDF. */
  Eigen::VectorXd _cdf;
};

}

#endif // __SBB_DIRECTMARGINALSAMPLER_HH__

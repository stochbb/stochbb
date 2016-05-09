#ifndef __SBB_EXACTSAMPLER_HH__
#define __SBB_EXACTSAMPLER_HH__

#include <Eigen/Eigen>
#include "api.hh"
#include <map>

namespace stochbb {

/** Implements an exact sampler for a system of random varaibles.
 * Proper sampling from a system of potentially dependent random variables is not trivial. First,
 * it must be ensured that every atomic random variable, a random variable that does not
 * depend on others, is sampled only once. Moreover, the samples of derived random variables need
 * to be cached for performance. Upon construction, this class first determines which and in
 * which order the random variables must be sampled. By calling @c sample, a sample is drawn from
 * the system.
 * @ingroup internal */
class ExactSamplerObj : public Object
{
public:
  /** Constructor for the given random variable. */
  ExactSamplerObj(const Var &X) throw (TypeError);
  /** Constructor for the given random variables. */
  ExactSamplerObj(const Var &X1, const Var &X2) throw (TypeError);
  /** Constructor for the given random variables. */
  ExactSamplerObj(const Var &X1, const Var &X2, const Var &X3) throw (TypeError);
  /** Constructor for the given random variables. */
  ExactSamplerObj(const std::vector<Var> &variables) throw (TypeError);

  virtual void mark();

  /** Samples the random variables given to the constructor.
   * The numer of rows of @c out specifies the number of samples, while the number of
   * columns of must match the number of variables passed to the constructor. */
  void sample(Eigen::Ref<Eigen::MatrixXd> out);

protected:
  /** Assembles the sampling queue recusrively. */
  void _add_to_queue(VarObj *var) throw (TypeError);

protected:
  /** The variables selected for output. */
  std::vector<VarObj *> _outvars;
  /** The queue of random variables to sample. */
  std::vector<VarObj *> _queue;
  /** A table mapping a random variable to its index in the queue. */
  std::map<VarObj *, size_t> _varmap;
};

}

#endif // __SBB_EXACTSAMPLER_HH__

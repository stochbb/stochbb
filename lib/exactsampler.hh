#ifndef __SBB_EXACTSAMPLER_HH__
#define __SBB_EXACTSAMPLER_HH__

#include <Eigen/Eigen>
#include "api.hh"
#include <map>

namespace sbb {

/** Implements an exact sampler for a system of random varaibles.
 * Proper sampling from a system of potentially dependent random variables is not trivial. First,
 * it must be ensured that every atomic random variable, a random variable that does not
 * depend on others, is sampled only once. Moreover, the samples of derived random variables need
 * to be cached for performance. Upon construction, this class first determines which and in
 * which order the random variables must be sampled. By calling @c sample, a sample is drawn from
 * the system. */
class ExactSamplerObj : public Object
{
public:
  /** Constructor from the given random variables. */
  ExactSamplerObj(const std::vector<Var> &variables);

  virtual void mark();

  /** Samples the random variables given to the constructor.
   * The numer of rows of @c out specifies the number of samples, while the number of
   * columns of must match the number of variables passed to the constructor. */
  void sample(Eigen::MatrixXd &out);

protected:
  /** Internal used fucntion type to sample a specific variable type. */
  typedef void (*sampler_f)(ExactSamplerObj *, VarObj *, Eigen::MatrixXd &out);
  /** Assembles the sampling queue recusrively. */
  void _add_to_queue(VarObj *var);
  /** Chooses the correct sampler for the given variable. */
  sampler_f _choose_sampler(VarObj *var);
  /** Samples from an atomic random variable. */
  static void _sample_atomic(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out);
  /** Samples from a chain of random variables. */
  static void _sample_chain(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out);
  /** Samples from a maximum of random variables. */
  static void _sample_minimum(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out);
  /** Samples from a minimum of random variables. */
  static void _sample_maximum(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out);
  /** Samples from a mixture of random variables. */
  static void _sample_mixture(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out);
  /** Samples from a compound-normal random variable. */
  static void _sample_comp_normal(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out);
  /** Samples from a compound-gamma random variable. */
  static void _sample_comp_gamma(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out);

protected:
  /** The variables selected for output. */
  std::vector<VarObj *> _outvars;
  /** The queue of random variables to sample. */
  std::vector<VarObj *> _queue;
  /** A table mapping a random variable to its index in the queue. */
  std::map<VarObj *, size_t> _varmap;
  /** The vector of samplers for each variable in the queue. */
  std::vector<sampler_f> _sampler;
};

}

#endif // __SBB_EXACTSAMPLER_HH__

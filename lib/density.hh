/** @defgroup density Densities
 * @ingroup internal */

#ifndef __SBB_DENSITY_HH__
#define __SBB_DENSITY_HH__

#include "object.hh"
#include "api.hh"
#include "exception.hh"
#include <Eigen/Eigen>

namespace stochbb {

/** Base class of all densities.
 * @ingroup density */
class DensityObj: public Object
{
protected:
  /** Hidden constructor. */
  DensityObj();

public:
  /** Destructor. */
  virtual ~DensityObj();
  void mark();

  /** Evaluates the density on a regular grid in \f$[Tmin, Tmax)\f$ and
   * stores it into the given output vector. The number of grid points is determined
   * by the length of the output vector.
   * @param Tmin Specifies the lower bound of the evaluation interval.
   * @param Tmax Specifies the upper bound of the evaluation interval.
   * @param out A preallocated vector that will hold the density evaluation. The length of the
   *   vector specifies the number of evaluation points. */
  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const = 0;

  /** Evaluates the probability function (CDF) at a regular grid on the interval \f$[Tmin, Tmax)\f$
   * and stores it into the given output vector. The number of grid points is determined
   * by the length of the output vector.
   * @param Tmin Specifies the lower bound of the evaluation interval.
   * @param Tmax Specifies the upper bound of the evaluation interval.
   * @param out A preallocated vector that will hold the CDF evaluation. The length of the
   *   vector specifies the number of evaluation points. */
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const = 0;

  /** Retruns a new density instance being the affine transformed of this density. That is
   * \f$f(x) \mapsto f(a\,x+b)\f$.
   * @param scale The scale parameter \f$a\f$.
   * @param shift The shift parameter \f$b\f$. */
  virtual Density affine(double scale, double shift) const = 0;

  /** Estimates the range on which the PDF and CDFs need to be evaluated based on the
   * \f$\alpha\f$-quantile of atomic distributions. Please note that this value
   * is not an estimate of \f$\alpha\f$-quantiles of the distribution. It reflects the interval on
   * which all distributions, this distribution depends on, need to be evaluated for the evaluation
   * of this distribution.
   * @param alpha Specifies the quantiles used for the estimation of the evaluation range.
   * @param a The lower bound estimate for the evaluation range.
   * @param b The upper bound estimate for the evaluation range. */
  virtual void rangeEst(double alpha, double &a, double &b) const = 0;

  /** Comparison operator between densities. This implementation compares only by type.
   * The specialization needs to compare also densities within types. */
  virtual int compare(const DensityObj &other) const;

  /** Prints a textual representation of the density object into the given @c stream. */
  virtual void print(std::ostream &stream) const;
};


/** Represents a specific "instantiation" of a distribution as a @c DensityObj for an atomic random
 * variable.
 * @ingroup density */
class AtomicDensityObj: public DensityObj
{
public:
  /** Constructs a generic atomic density with the given distribution and parameters.
   * @param dist Specifies the distribution of the random variable.
   * @param params Specifies the parameters for the random variable. */
  AtomicDensityObj(const Distribution &dist, Eigen::Ref<Eigen::VectorXd> params);
  /** Destructor. */
  virtual ~AtomicDensityObj();
  void mark();

  /** Returns a reference to the distribution object. */
  Distribution distribution() const;
  /** Returns the number of parameters of the density. */
  size_t nParams() const;
  /** Returns a specific parameter. */
  double parameter(size_t i) const;

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;

  Density affine(double scale, double shift) const;

  void rangeEst(double alpha, double &a, double &b) const;

  /** Samples from the atomic distribution. */
  void sample(Eigen::Ref<Eigen::VectorXd> out) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

protected:
  /** Holds a reference to the distribution instance. */
  DistributionObj *_distribution;
  /** Holds the parameter vector. */
  Eigen::VectorXd _params;
};


}

#endif // __SBB_DENSITY_HH__

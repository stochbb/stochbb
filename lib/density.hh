#ifndef __SBB_DENSITY_HH__
#define __SBB_DENSITY_HH__

#include "object.hh"
#include "api.hh"
#include <Eigen/Eigen>

namespace stochbb {


/** Base class of all densities. */
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


/** Implements the base class of all densities of atomic random variables.
 * As atomic variables do not depend on other random variables, it is possible to sample them
 * directly. Hence their densities implement a @c sample method. */
class AtomicDensityObj: public DensityObj
{
protected:
  /** Hidden constructor. */
  AtomicDensityObj();

public:
  /** Destructor. */
  virtual ~AtomicDensityObj();

  void mark();

  /** Samples from the atomic density. Only atomic densities, that is densities which do not depend
   * on other densities, can be sampled directly. For all other densities consider the
   * @c ExactSampler or @c MarginalSampler classes. */
  virtual void sample(Eigen::Ref<Eigen::VectorXd> out) const = 0;
};


/** Implements the delta distribution. */
class DeltaDensityObj: public AtomicDensityObj
{
public:
  /** Constructor. */
  DeltaDensityObj(double delay);
  /** Destructor. */
  virtual ~DeltaDensityObj();
  void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out) const;
  Density affine(double scale, double shift) const;
  void rangeEst(double alpha, double &a, double &b) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

  /** Retruns the delay of the delta distribution. */
  inline double delay() const { return _delay; }

protected:
  /** Holds the center of the distribution. */
  double _delay;
};


/** Implements the uniform distribution on the interval \f$[a,b]\f$. */
class UniformDensityObj: public AtomicDensityObj
{
public:
  /** Constructs a uniform distribution on the interval \f$[a,b]\f$.
   * @param a Specifies the lower bound of the interval.
   * @param b Specifies the upper bound of the interval. */
  UniformDensityObj(double a, double b);
  /** Destructor. */
  virtual ~UniformDensityObj();
  void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out) const;
  Density affine(double scale, double shift) const;
  void rangeEst(double alpha, double &a, double &b) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

  /** Returns the lower-bound of the uniform distribution. */
  inline double a() const { return _a; }
  /** Returns the upper-bound of the uniform distribution. */
  inline double b() const { return _b; }

protected:
  /** The lower end of the interval. */
  double _a;
  /** The upper end of the interval. */
  double _b;
};


/** Implements the normal distribution.
 * That is
 * \f[
 *  \phi(x, \mu, \sigma) = \frac{\exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)}{\sqrt{2\pi}\sigma}\,.
 * \f]
 */
class NormalDensityObj: public AtomicDensityObj
{
public:
  /** Constructor with mean and standard deviation. */
  NormalDensityObj(double mean, double stddev);
  /** Destructor. */
  virtual ~NormalDensityObj();
  void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out) const;
  Density affine(double scale, double shift) const;
  void rangeEst(double alpha, double &a, double &b) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

  /** Returns the mean of the normal distribution. */
  inline double mu() const { return _mu; }
  /** Returns the standard deviation of the normal distribution. */
  inline double sigma() const { return _sigma; }

protected:
  /** The mean. */
  double _mu;
  /** The standard deviation. */
  double _sigma;
};


/** Implements the (shifted) Gamma distribution. That is
 * \f[
 *  \Gamma(x; k, \theta) = \frac{x^{k-1}e^{-\frac{x}{\theta}}}{\Gamma(k)\,\theta^k}
 * \f] */
class GammaDensityObj: public AtomicDensityObj
{
public:
  /** Constructor. */
  GammaDensityObj(double k, double theta, double shift=0);
  /** Destructor. */
  virtual ~GammaDensityObj();
  void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out) const;
  Density affine(double scale, double shift) const;
  void rangeEst(double alpha, double &a, double &b) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

  /** Returns the shape parameter of the gamma distribution. */
  inline double k() const { return _k; }
  /** Returns the scale parameter of the gamma distribution. */
  inline double theta() const { return _theta; }
  /** Returns the shift of the distribution. */
  inline double shift() const { return _shift; }

protected:
  /** The shape paramter. */
  double _k;
  /** The scale parameter. */
  double _theta;
  /** Shift of the affine transform. */
  double _shift;
};


/** Implements the (shifted) inverse Gamma distribution. That is
 * \f[
 *  \Gamma(x; \alpha, \beta) = \frac{\beta^\alpha x^{-\alpha-1}e^{-\frac{\beta}{a}}}{\Gamma(\alpha)}
 * \f] */
class InvGammaDensityObj: public AtomicDensityObj
{
public:
  /** Constructor. */
  InvGammaDensityObj(double alpha, double beta, double shift=0);
  /** Destructor. */
  virtual ~InvGammaDensityObj();
  void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out) const;
  Density affine(double scale, double shift) const;
  void rangeEst(double alpha, double &a, double &b) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

  /** Returns the shape parameter of the inverse gamma distribution. */
  inline double alpha() const { return _alpha; }
  /** Returns the scale parameter of the inverse gamma distribution. */
  inline double beta() const { return _beta; }
  /** Returns the shift of the distribution. */
  inline double shift() const { return _shift; }

protected:
  /** The shape paramter. */
  double _alpha;
  /** The scale parameter. */
  double _beta;
  /** Shift of the affine transform. */
  double _shift;
};


/** Implements the Weibull distribution.
 * The (shifted) Weilbull distribution is defined as
 * \f[
 *  f(x\ge 0;k,\lambda,\theta) = \frac{k}{\lambda}\left(\frac{x-\theta}{\lambda}\right)^{k-1}
 *    \exp\left[-\left(\frac{x-\theta}{\lambda}\right)^k\right]\,,
 * \f]
 * where \f$k\f$ is the shape parameter, \f$\lambda\f$ the scale parameter and \f$\theta\f$
 * the shift*/
class WeibullDensityObj: public AtomicDensityObj
{
public:
  /** Constructor with shape @c k, scale @c lambda and @c shift. */
  WeibullDensityObj(double k, double lambda, double shift=0);
  /** Destructor. */
  virtual ~WeibullDensityObj();
  void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out) const;
  Density affine(double scale, double shift) const;
  void rangeEst(double alpha, double &a, double &b) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

  /** Returns the shape parameter of the Weilbull distribution. */
  inline double k() const { return _k; }
  /** Returns the scale parameter of the Weibull distribution. */
  inline double lambda() const { return _lambda; }
  /** Returns the shift of the distribution. */
  inline double shift() const { return _shift; }

protected:
  /** The shape paramter. */
  double _k;
  /** The scale parameter. */
  double _lambda;
  /** Shift of the affine transform. */
  double _shift;
};

}

#endif // __SBB_DENSITY_HH__

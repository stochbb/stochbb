#ifndef __SSB_DISTRIBUTIONOBJ_HH__
#define __SSB_DISTRIBUTIONOBJ_HH__

#include "density.hh"
#include <Eigen/Eigen>


namespace stochbb {

/** Base class of all distributions, that are families of distribution functions that are
 * parametrized by a set of parameters. An probability (atomic) density that is assigned to a
 * random variable can then be considered as an instance of that family with a specific set of
 * parameters (@c GenericAtomicDensityObj).
 *
 * All densities have at least 2 parameters which are \f$[scale, shift, ...]\f$ such that a specific
 * density can be expressed as \f$f\left(\frac{x-shift}{scale}, ...\right)\f$. */
class DistributionObj: public Object
{
protected:
  /** Hidden constructor. */
  DistributionObj();

public:
  /** Destructor. */
  virtual ~DistributionObj();

  void mark();

  /** Returns the number of parameters. */
  virtual size_t nParams() const = 0;
  /** Evaluates the probability density function at the given point @c x with the specified parameters. */
  virtual double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Evaluates the probability function at the given point @c x with the specified parameters. */
  virtual double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Returns the quantiles (@c lower, @c upper) for the given probability @c p with the specified
   * parameters. */
  virtual void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Changes the given set of parameters such that the distribution is an affine transformed. */
  virtual void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const = 0;
  /** Changes the given set of parameter densities such that the distribution is an affine
   * transformed. */
  virtual void affine(double scale, double shift, std::vector<DensityObj *> &params) const = 0;
  /** Draws a sample from the distribution with the specified parameters. */
  virtual double sample(const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Distribution type order. */
  virtual int compare(const DistributionObj &other) const;
};


/** The family of uniform densities. */
class UniformDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  UniformDistributionObj();
  /** Destructor. */
  virtual ~UniformDistributionObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  void affine(double scale, double shift, std::vector<DensityObj *> &params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


/** The family of normal densities. */
class NormalDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  NormalDistributionObj();
  /** Destructor. */
  virtual ~NormalDistributionObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


/** The family of gamma densities. */
class GammaDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  GammaDistributionObj();
  /** Destructor. */
  virtual ~GammaDistributionObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


/** The family of inverse gamma densities. */
class InvGammaDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  InvGammaDistributionObj();
  /** Destructor. */
  virtual ~InvGammaDistributionObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


/** The family of Weibull densities. */
class WeibullDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  WeibullDistributionObj();
  /** Destructor. */
  virtual ~WeibullDistributionObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


/** Represents a specific "instantiation" of a distribution as a @c DensityObj of atomic random
 * variables. */
class GenericAtomicDensityObj: public AtomicDensityObj
{
public:
  /** Constructs a generic atomic density with the given distribution and parameters. */
  GenericAtomicDensityObj(DistributionObj *dist, Eigen::Ref<Eigen::VectorXd> params);
  /** Destructor. */
  virtual ~GenericAtomicDensityObj();
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

  void sample(Eigen::Ref<Eigen::VectorXd> out) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

protected:
  /** Holds a reference to the distribution instance. */
  DistributionObj *_distribution;
  /** Holds the parameter vector. */
  Eigen::VectorXd _params;
};


/** Implements the generic density of compound random variables. That is a random variable
 * that is distributed with an atomic density where the parameters are random variables too. */
class GenericCompoundDensityObj: public DensityObj
{
public:
  /** Constructs a compound density with the given distribution family and parameter distributions. */
  GenericCompoundDensityObj(DistributionObj *dist, const std::vector<DensityObj *> &params);
  /** Destructor. */
  virtual ~GenericCompoundDensityObj();
  void mark();

  /** Returns a reference to the distribution object. */
  Distribution distribution() const;
  /** Returns the number of parameters of the density. */
  size_t nParams() const;
  /** Returns a specific parameter density. */
  Density parameter(size_t i) const;

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;

  Density affine(double scale, double shift) const;

  void rangeEst(double alpha, double &a, double &b) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

protected:
  /** Helper function to iterate over several indices.
   * @param i Specifies the running index.
   * @param N Specifies the number of dimensions per index.
   * @param idxs On exit, contains the vector of indices. */
  void _to_param_indices(size_t i, size_t N, std::vector<size_t> &idxs) const;

protected:
  /** Holds a reference to the distribution instance. */
  DistributionObj *_distribution;
  /** Holds references to the parameter densities. */
  std::vector<DensityObj *> _parameters;
};

}

#endif // __SSB_DISTRIBUTIONOBJ_HH__

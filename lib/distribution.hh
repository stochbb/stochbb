#ifndef __SSB_DISTRIBUTIONOBJ_HH__
#define __SSB_DISTRIBUTIONOBJ_HH__

#include "density.hh"
#include <Eigen/Eigen>


namespace stochbb {

/** Base class of all distributions, that are families of distribution functions that are
 * parametrized by a set of parameters. An probability (atomic) density that is assigned to a
 * random variable can then be considered as an instance of that family with a specific set of
 * parameters (@c GenericAtomicDensityObj). */
class DistributionObj: public Object
{
protected:
  /** Hidden constructor. */
  DistributionObj();

public:
  /** Destructor. */
  virtual ~DistributionObj();

  void mark();

  /** Return the number of parameters. */
  virtual size_t nParams() const = 0;
  /** Evaluates the probability density function at the given point @c x with the specified parameters. */
  virtual double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Evaluates the probability function at the given point @c x with the specified parameters. */
  virtual double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Returns the quantiles (@c lower, @c upper) for the given probability @c p with the specified
   * parameters. */
  virtual void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Changes the given set of parameters being the affine transformed of this distribution. */
  virtual void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const = 0;
  /** Changes the given set of parameter random variables being the affine transformed of this
   * distribution. */
  virtual void affine(double scale, double shift, std::vector<DensityObj *> &params) const = 0;
  /** Draws a sample from the distribution with the specified parameters. */
  virtual double sample(const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Distribution type order. */
  virtual int compare(const DistributionObj &other) const;
};


class GenericAtomicDensityObj: public AtomicDensityObj
{
public:
  GenericAtomicDensityObj(DistributionObj *dist, Eigen::Ref<Eigen::VectorXd> params);
  virtual ~GenericAtomicDensityObj();
  void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;

  Density affine(double scale, double shift) const;

  void rangeEst(double alpha, double &a, double &b) const;

  void sample(Eigen::Ref<Eigen::VectorXd> out) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

protected:
  DistributionObj *_distribution;
  Eigen::VectorXd _params;
};


class UniformDistributionObj: public DistributionObj
{
public:
  UniformDistributionObj();
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


class NormalDistributionObj: public DistributionObj
{
public:
  NormalDistributionObj();
  virtual ~NormalDistributionObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


class GammaDistributionObj: public DistributionObj
{
public:
  GammaDistributionObj();
  virtual ~GammaDistributionObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


class InvGammaDistributionObj: public DistributionObj
{
public:
  InvGammaDistributionObj();
  virtual ~InvGammaDistributionObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


class WeibullDistributionObj: public DistributionObj
{
public:
  WeibullDistributionObj();
  virtual ~WeibullDistributionObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


class GenericCompoundDensityObj: public DensityObj {
public:
  GenericCompoundDensityObj(DistributionObj *dist, const std::vector<DensityObj *> &params);
  virtual ~GenericCompoundDensityObj();
  void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;

  Density affine(double scale, double shift) const;

  void rangeEst(double alpha, double &a, double &b) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

protected:
  void _to_param_indices(size_t i, size_t N, std::vector<size_t> &idxs) const;

protected:
  DistributionObj *_distribution;
  std::vector<DensityObj *> _parameters;
};

}

#endif // __SSB_DISTRIBUTIONOBJ_HH__

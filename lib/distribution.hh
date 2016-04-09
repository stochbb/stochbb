#ifndef __SSB_DISTRIBUTIONOBJ_HH__
#define __SSB_DISTRIBUTIONOBJ_HH__

#include "density.hh"
#include <Eigen/Eigen>


namespace stochbb {

/** Base class of all distributions. */
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
  /** Evaluates the density function at the given point @c x with the specified parameters. */
  virtual double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Evaluates the probability function at the given point @c x with the specified parameters. */
  virtual double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Returns the quantile for the give probability @c p with the specified parameters. */
  virtual void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Changes the given set of parameters being the affine transformed of this distribution. */
  virtual void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const = 0;
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


class UniformDistObj: public DistributionObj
{
public:
  UniformDistObj();
  virtual ~UniformDistObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


class NormalDistObj: public DistributionObj
{
public:
  NormalDistObj();
  virtual ~NormalDistObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


class GammaDistObj: public DistributionObj
{
public:
  GammaDistObj();
  virtual ~GammaDistObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


class InvGammaDistObj: public DistributionObj
{
public:
  InvGammaDistObj();
  virtual ~InvGammaDistObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};


class WeibullDistObj: public DistributionObj
{
public:
  WeibullDistObj();
  virtual ~WeibullDistObj();
  void mark();

  size_t nParams() const;
  double pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;
};

}

#endif // __SSB_DISTRIBUTIONOBJ_HH__

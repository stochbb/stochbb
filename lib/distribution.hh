#ifndef __SSB_DISTRIBUTIONOBJ_HH__
#define __SSB_DISTRIBUTIONOBJ_HH__

#include "density.hh"
#include <Eigen/Eigen>


namespace stochbb {

/** Base class of all distributions. Distributions are families of density functions that are
 * parametrized by a set of parameters. A (atomic) probability density, that is assigned to a
 * random variable can then be considered as an instance of that family with a specific set of
 * parameters (@c GenericAtomicDensityObj).
 *
 * Almost all densities (except delta) have at least 2 parameters which are \f$[scale, shift, ...]\f$
 * such that a specific density can be expressed as \f$f\left(\frac{x-shift}{scale}, ...\right)\f$.
 * @see AtomicDensity @see Compound
 * @ingroup density */
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
  /** Evaluates the probability density function with the specified parameters on a regular grid
   * in \f$[T_{min}, T_{max})\f$ and stores the result into @c out. */
  virtual void pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
                   const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Evaluates the probability function with the specified parameters on a regular grid
   * in \f$[T_{min}, T_{max})\f$ and stores the result into @c out. */
  virtual void cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
                   const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Returns the quantiles (@c lower, @c upper) for the given probability @c p with the specified
   * parameters. */
  virtual void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Changes the given set of parameters such that the distribution is an affine transformed. */
  virtual void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const = 0;
  /** Changes the given set of parameter densities such that the distribution is an affine
   * transformed. */
  virtual void affine(double scale, double shift, std::vector<Density> &params) const = 0;
  /** Draws some samples from the distribution with the specified parameters. */
  virtual void sample(Eigen::Ref<Eigen::VectorXd> out,
                      const Eigen::Ref<const Eigen::VectorXd> params) const = 0;
  /** Distribution type order. */
  virtual int compare(const DistributionObj &other) const;
};


/** The family of delta distributions.
 * @ingroup density */
class DeltaDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  DeltaDistributionObj();
  /** Destructor. */
  virtual ~DeltaDistributionObj();
  void mark();

  size_t nParams() const;
  void pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  void affine(double scale, double shift, std::vector<Density> &params) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void print(std::ostream &stream) const;
};


/** The family of uniform distributions.
 * @ingroup density */
class UniformDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  UniformDistributionObj();
  /** Destructor. */
  virtual ~UniformDistributionObj();
  void mark();

  size_t nParams() const;
  void pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  void affine(double scale, double shift, std::vector<Density> &params) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void print(std::ostream &stream) const;
};


/** The family of normal distributions.
 * @ingroup density */
class NormalDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  NormalDistributionObj();
  /** Destructor. */
  virtual ~NormalDistributionObj();
  void mark();

  size_t nParams() const;
  void pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  void affine(double scale, double shift, std::vector<Density> &params) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void print(std::ostream &stream) const;
};


/** The family of gamma distributions.
 * @ingroup density */
class GammaDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  GammaDistributionObj();
  /** Destructor. */
  virtual ~GammaDistributionObj();
  void mark();

  size_t nParams() const;
  void pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  void affine(double scale, double shift, std::vector<Density> &params) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void print(std::ostream &stream) const;
};


/** The family of inverse gamma distributions.
 * @ingroup density */
class InvGammaDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  InvGammaDistributionObj();
  /** Destructor. */
  virtual ~InvGammaDistributionObj();
  void mark();

  size_t nParams() const;
  void pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  void affine(double scale, double shift, std::vector<Density> &params) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void print(std::ostream &stream) const;
};


/** The family of Weibull distributions.
 * @ingroup density */
class WeibullDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  WeibullDistributionObj();
  /** Destructor. */
  virtual ~WeibullDistributionObj();
  void mark();

  size_t nParams() const;
  void pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  void affine(double scale, double shift, std::vector<Density> &params) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void print(std::ostream &stream) const;
};


/** The family of (non-standardized) Student's t distributions.
 * @ingroup density */
class StudtDistributionObj: public DistributionObj
{
public:
  /** Constructor. */
  StudtDistributionObj();
  virtual ~StudtDistributionObj();

  size_t nParams() const;
  void pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  void affine(double scale, double shift, std::vector<Density> &params) const;
  void sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const;
  void print(std::ostream &stream) const;
};

}

#endif // __SSB_DISTRIBUTIONOBJ_HH__

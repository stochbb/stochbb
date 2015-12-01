#ifndef __SBB_DENSITY_HH__
#define __SBB_DENSITY_HH__

#include "object.hh"
#include <Eigen/Eigen>

namespace sbb {


/** Base class of all densities.
 * @ingroup internal */
class DensityObj: public Object
{
protected:
  /** Hidden constructor. */
  DensityObj();

public:
  /** Destructor. */
  virtual ~DensityObj();
  virtual void mark();

  /** Evaluates the density at a regular grid on the interval \f$[Tmin, Tmax)\f$ and
   * stores it into the given output vector. The number of grid points is determined
   * by the length of the output vector. */
  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const = 0;
  /** Evaluates the probability function (CDF) at a regular grid on the interval \f$[Tmin, Tmax)\f$
   * and stores it into the given output vector. The number of grid points is determined
   * by the length of the output vector. */
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const = 0;
  /** Samples from the distribution. */
  virtual void sample(Eigen::VectorXd &out) const = 0;
};


/** Implements the delta distribution.
 * @ingroup internal */
class DeltaDensityObj: public DensityObj
{
public:
  /** Constructor. */
  DeltaDensityObj(double delay);
  /** Destructor. */
  virtual ~DeltaDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

protected:
  /** Holds the center of the distribution. */
  double _delay;
};


/** Implements the uniform distribution on the interval \f$[a,b]\f$.
 * @ingroup internal */
class UniformDensityObj: public DensityObj
{
public:
  /** Constructor. */
  UniformDensityObj(double a, double b);
  /** Destructor. */
  virtual ~UniformDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

protected:
  /** The lower end of the interval. */
  double _a;
  /** The upper end of the interval. */
  double _b;
};


/** Implements the normal distribution.
 * @ingroup internal */
class NormalDensityObj: public DensityObj
{
public:
  /** Constructor with mean and standard deviation. */
  NormalDensityObj(double mean, double stddev);
  /** Destructor. */
  virtual ~NormalDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

protected:
  /** The mean. */
  double _mu;
  /** The standard deviation. */
  double _sigma;
};


/** Implements the Gamma distribution.
 * @ingroup internal */
class GammaDensityObj: public DensityObj
{
public:
  /** Constructor. */
  GammaDensityObj(double k, double theta);
  /** Destructor. */
  virtual ~GammaDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

protected:
  /** The shape paramter. */
  double _k;
  /** The scale parameter. */
  double _theta;
};

}

#endif // __SBB_DENSITY_HH__

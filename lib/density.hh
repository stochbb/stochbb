#ifndef __SBB_DENSITY_HH__
#define __SBB_DENSITY_HH__

#include "object.hh"
#include "api.hh"
#include <Eigen/Eigen>

namespace sbb {


/** Base class of all densities. */
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

  /** Retruns a new density instance being the affine transformed of this density. */
  virtual Density affine(double scale, double shift) const = 0;

  /** Comparison operator between densities. This implementation compares only by type.
   * The specialization needs to compare also densities within types. */
  virtual int compare(const DensityObj &other) const;

  /** Prints a textual representation of the density object. */
  virtual void print(std::ostream &stream) const;
};

/** Prints a textual representation of the density.
 * @ingroup internal */
inline std::ostream &operator<<(std::ostream &stream, const DensityObj &density) {
  density.print(stream); return stream;
}


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

  virtual void mark();

  /** Samples from the density. */
  virtual void sample(Eigen::VectorXd &out) const = 0;
};


/** Implements the delta distribution. */
class DeltaDensityObj: public AtomicDensityObj
{
public:
  /** Constructor. */
  DeltaDensityObj(double delay, double scale=1, double shift=0);
  /** Destructor. */
  virtual ~DeltaDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;
  virtual Density affine(double scale, double shift) const;

  /** Compares densities. */
  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

  /** Retruns the delay of the delta distribution. */
  inline double delay() const { return _delay; }

protected:
  /** Holds the center of the distribution. */
  double _delay;
  /** Holds the constants of an affine tranformation of this density. */
  double _scale;
  /** Holds the constants of an affine tranformation of this density. */
  double _shift;
};


/** Implements the uniform distribution on the interval \f$[a,b]\f$. */
class UniformDensityObj: public AtomicDensityObj
{
public:
  /** Constructor. */
  UniformDensityObj(double a, double b, double scale=1, double shift=0);
  /** Destructor. */
  virtual ~UniformDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;
  virtual Density affine(double scale, double shift) const;

  /** Compares densities. */
  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

  /** Returns the lower-bound of the uniform distribution. */
  inline double a() const { return _a; }
  /** Returns the upper-bound of the uniform distribution. */
  inline double b() const { return _b; }

protected:
  /** The lower end of the interval. */
  double _a;
  /** The upper end of the interval. */
  double _b;
  /** Scaling of the affine transformation. */
  double _scale;
  /** Shift of the affine transformation. */
  double _shift;
};


/** Implements the normal distribution. */
class NormalDensityObj: public AtomicDensityObj
{
public:
  /** Constructor with mean and standard deviation. */
  NormalDensityObj(double mean, double stddev, double scale=1, double shift=0);
  /** Destructor. */
  virtual ~NormalDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;
  virtual Density affine(double scale, double shift) const;

  /** Compares densities. */
  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

  /** Returns the mean of the normal distribution. */
  inline double mu() const { return _mu; }
  /** Returns the standard deviation of the normal distribution. */
  inline double sigma() const { return _sigma; }

protected:
  /** The mean. */
  double _mu;
  /** The standard deviation. */
  double _sigma;
  /** Scale of the affine transformation. */
  double _scale;
  /** Shift of the affine transformation. */
  double _shift;
};


/** Implements the Gamma distribution. */
class GammaDensityObj: public AtomicDensityObj
{
public:
  /** Constructor. */
  GammaDensityObj(double k, double theta, double scale=1, double shift=0);
  /** Destructor. */
  virtual ~GammaDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;
  virtual Density affine(double scale, double shift) const;

  /** Compares densities. */
  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

  /** Returns the shape parameter of the gamma distribution. */
  inline double k() const { return _k; }
  /** Returns the scale parameter of the gamma distribution. */
  inline double theta() const { return _theta; }

protected:
  /** The shape paramter. */
  double _k;
  /** The scale parameter. */
  double _theta;
  /** Scale of the affine transform. */
  double _scale;
  /** Shift of the affine transform. */
  double _shift;
};

}

#endif // __SBB_DENSITY_HH__

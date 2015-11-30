#ifndef __SBB_DENSITY_HH__
#define __SBB_DENSITY_HH__

#include "object.hh"
#include <Eigen/Eigen>

namespace sbb {


class DensityObj: public Object
{
protected:
  DensityObj();

public:
  virtual ~DensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const = 0;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const = 0;
  virtual void sample(Eigen::VectorXd &out) const = 0;
};


class DeltaDensityObj: public DensityObj
{
public:
  DeltaDensityObj(double delay);
  virtual ~DeltaDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

protected:
  double _delay;
};


class UniformDensityObj: public DensityObj
{
public:
  UniformDensityObj(double a, double b);
  virtual ~UniformDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

protected:
  double _a;
  double _b;
};


class NormalDensityObj: public DensityObj
{
public:
  NormalDensityObj(double mean, double stddev);
  virtual ~NormalDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

protected:
  double _mu;
  double _sigma;
};


class GammaDensityObj: public DensityObj
{
public:
  GammaDensityObj(double k, double theta);
  virtual ~GammaDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

protected:
  double _k;
  double _theta;
};

}

#endif // __SBB_DENSITY_HH__

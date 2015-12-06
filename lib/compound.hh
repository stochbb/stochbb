#ifndef __SBB_COMPOUND_HH__
#define __SBB_COMPOUND_HH__

#include "api.hh"
#include "randomvariable.hh"

namespace sbb {

/** Base class of all compound random variables. */
class CompoundObj : public DerivedVarObj
{
protected:
  /** Hidden constructor.
   * @param vars Specifies the variables, the compound depends on.
   * @param name Specifies the optional name for the random variable. */
  CompoundObj(const std::vector<Var> &vars, const std::string &name);

public:
  /** Constructs a compound-normal random variable with the specified random variables
   * for the mean and standard deviation. */
  static Compound norm(const Var &mu, const Var &sigma, const std::string &name="");
  /** Constructs a compound-gamma random variable with the specified random variables
   * for the shape and scale parameters. */
  static Compound gamma(const Var &k, const Var &theta, const std::string &name="");

public:
  virtual void mark();
};


/** The density object for a compound-normal distribution.
 * This class determines the distribution of the random variable numerically. Depending on the
 * number of gird points at which the parameter distributions are evaluated, this may take a
 * considerable time. */
class NormalCompoundDensityObj: public DensityObj
{
public:
  /** Constructs the density from the given parameter random variables. */
  NormalCompoundDensityObj(const Var &mu, const Var &sigma);

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

protected:
  /** A reference to the distribution of the mean parameter. */
  DensityObj *_mu;
  /** A reference to the distribution of the standard deviation parameter. */
  DensityObj *_sigma;
};


/** Implements a compound-normal distribution. */
class NormalCompoundObj: public CompoundObj
{
public:
  /** Constructs a compound-normal distribution from the given parameter random variables. */
  NormalCompoundObj(const Var &mu, const Var &sigma, const std::string &name="");

  void mark();

  virtual Density density();

protected:
  /** A reference to the density. */
  NormalCompoundDensityObj *_density;
};


/** The density object for a compound-gamma distribution.
 * This class determines the distribution of the random variable numerically. Depending on the
 * number of gird points at which the parameter distributions are evaluated, this may take a
 * considerable time. */
class GammaCompoundDensityObj: public DensityObj
{
public:
  /** Constructs the density from the given parameter distributions. */
  GammaCompoundDensityObj(const Var &k, const Var &theta);

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

protected:
  /** The distribution of the shape parameter. */
  DensityObj *_k;
  /** The distribution of the scale parameter. */
  DensityObj *_theta;
};


/** Implements a compound-gamma random variable. */
class GammaCompoundObj: public CompoundObj
{
public:
  /** Constructs the compound-gamma random variable from the given parameter random variables. */
  GammaCompoundObj(const Var &k, const Var &theta, const std::string &name="");

  void mark();

  virtual Density density();

protected:
  /** A reference to the distribution. */
  GammaCompoundDensityObj *_density;
};

}

#endif // COMPOUND_HH

#ifndef __SBB_COMPOUND_HH__
#define __SBB_COMPOUND_HH__

#include "api.hh"
#include "randomvariable.hh"

namespace stochbb {

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
protected:
  /** Constructs the density from the given parameter random variables.
   * @param mu Specifies the mean random variable.
   * @param sigma Specifies the standard deviation random variable.
   * @param scale Specifies the optional scaling, default @c scale=1.
   * @param shift Specifies the optional shift, default @c shift=0. */
  NormalCompoundDensityObj(DensityObj *mu, DensityObj *sigma, double scale=1, double shift=0);

public:
  /** Constructs the density from the given parameter random variables.
   * @param mu Specifies the mean random variable.
   * @param sigma Specifies the standard deviation random variable.
   * @param scale Specifies the optional scaling, default @c scale=1.
   * @param shift Specifies the optional shift, default @c shift=0. */
  NormalCompoundDensityObj(const Var &mu, const Var &sigma, double scale=1, double shift=0);

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;

protected:
  /** A reference to the distribution of the mean parameter. */
  DensityObj *_mu;
  /** A reference to the distribution of the standard deviation parameter. */
  DensityObj *_sigma;
  /** Scale of the affine transform. */
  double _scale;
  /** Shift of the affine transform. */
  double _shift;
};


/** Implements a compound-normal distribution.
 * That is a normal-distributed random variable where the mean and standard deviation are random
 * variables too. */
class NormalCompoundObj: public CompoundObj
{
public:
  /** Constructs a compound-normal distribution from the given parameter random variables.
   * @param mu Specifies the mean random variable.
   * @param sigma Specifies the standard deviation random variable.
   * @param name The optional variable name. */
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
protected:
  /** Constructs the density from the given parameter distributions.
   * @param k Specifies the shape parameter random variable.
   * @param theta Specifies the scale parameter random variable.
   * @param scale Specifies the optional scaling, default @c scale=1.
   * @param shift Specifies the optional shift, default @c shift=0. */
  GammaCompoundDensityObj(DensityObj *k, DensityObj *theta, double scale=1, double shift=0);

public:
  /** Constructs the density from the given parameter distributions.
   * @param k Specifies the shape parameter random variable.
   * @param theta Specifies the scale parameter random variable.
   * @param scale Specifies the optional scaling, default @c scale=1.
   * @param shift Specifies the optional shift, default @c shift=0. */
  GammaCompoundDensityObj(const Var &k, const Var &theta, double scale=1, double shift=0);

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;

protected:
  /** The distribution of the shape parameter. */
  DensityObj *_k;
  /** The distribution of the scale parameter. */
  DensityObj *_theta;
  /** Scale of the affine transform. */
  double _scale;
  /** Shift of the affine transform. */
  double _shift;
};


/** Implements a compound-gamma random variable.
 * That is a gamma-distributed random variable where the shape and scale parameters are random
 * variables too. */
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

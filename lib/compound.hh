#ifndef __SBB_COMPOUND_HH__
#define __SBB_COMPOUND_HH__

#include "api.hh"
#include "randomvariable.hh"
#include <Eigen/Eigen>


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
  virtual void mark();
};


/** The density object for a compound-normal distribution.
 * This class determines the distribution of the random variable numerically. That is
 * \f[
 *  f(x) = \int d\mu \int d\sigma\,\phi(x;\mu, \sigma)\,g(\mu)\,h(\sigma)\,,
 * \f]
 * where \f$g(\mu)\f$ and \f$h(\sigma)\f$ are the densities of the mean and std. dev. parameters,
 * repsectively. Numerically, the intergals are approximated by a finite sum, where the integration
 * interval is determined by the @c DensityObj::rangeEst methods of the parameter densities using
 * \f$\alpha=0.01\f$. */
class NormalCompoundDensityObj: public DensityObj
{
protected:
  /** Constructs the compound density from the given parameter random variables.
   * @param mu Specifies the mean random variable.
   * @param sigma Specifies the standard deviation random variable.*/
  NormalCompoundDensityObj(DensityObj *mu, DensityObj *sigma);

public:
  /** Constructs the compound density from the given random variables.
   * @param mu Specifies the mean random variable.
   * @param sigma Specifies the standard deviation random variable.
   * @throws AssumptionError If the parameter variables are not mutually independent. */
  NormalCompoundDensityObj(const Var &mu, const Var &sigma) throw (AssumptionError);
  /** Destructor. */
  virtual ~NormalCompoundDensityObj();

  void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  Density affine(double scale, double shift) const;
  void rangeEst(double alpha, double &a, double &b) const;

protected:
  /** Prepares the numerical itergration on construction. If one of the parameter densities is a
   * delta density, the integration over this parameter is performed "analytically". */
  void _init_int();

protected:
  /** A reference to the distribution of the mean parameter. */
  DensityObj *_mu;
  /** A reference to the distribution of the standard deviation parameter. */
  DensityObj *_sigma;
  /** Lower bound of the integration interval over \f$\mu\f$. */
  double _muMin;
  /** Integration stepsize over \f$\mu\f$. */
  double _ddMu;
  /** The density of \f$\mu\f$ evaluated on a regular grid from @c _muMin with step size @c _ddMu. */
  Eigen::VectorXd _dmu;
  /** Lower bound of the integration interval over \f$\sigma\f$. */
  double _sigMin;
  /** Integration stepsize over \f$\sigma\f$. */
  double _ddSig;
  /** The density of \f$\sigma\f$ evaluated on a regular grid from @c _sigMin with step size @c _ddSig. */
  Eigen::VectorXd _dsigma;
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
 * This class determines the distribution of the random variable numerically by integrating over
 * the densities of the specified parameter random variables.
 */
class GammaCompoundDensityObj: public DensityObj
{
protected:
  /** Constructs the density from the given parameter distributions.
   * @param k Specifies the shape parameter random variable.
   * @param theta Specifies the scale parameter random variable.
   * @param shift Specifies the optional shift, default @c shift=0. */
  GammaCompoundDensityObj(DensityObj *k, DensityObj *theta, double shift=0);

public:
  /** Constructs the density from the given parameter distributions.
   * @param k Specifies the shape parameter random variable.
   * @param theta Specifies the scale parameter random variable.
   * @param shift Specifies the optional shift, default @c shift=0. */
  GammaCompoundDensityObj(const Var &k, const Var &theta, double shift=0) throw (AssumptionError);
  virtual ~GammaCompoundDensityObj();

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;
  virtual void rangeEst(double alpha, double &a, double &b) const;

protected:
  /** Prepares the numerical itergration on construction. If one of the parameter densities is a
   * delta density, the integration over this parameter is performed "analytically". */
  void _init_int();

protected:
  /** The distribution of the shape parameter. */
  DensityObj *_k;
  /** The distribution of the scale parameter. */
  DensityObj *_theta;
  /** Shift of the affine transform. */
  double _shift;
  /** Lower bound of the integration interval over \f$k\f$. */
  double _kMin;
  /** Integration stepsize over \f$k\f$. */
  double _ddK;
  /** The density of \f$k\f$ evaluated on a regular grid from @c _kMin with step size @c _ddK. */
  Eigen::VectorXd _dk;
  /** Lower bound of the integration interval over \f$\theta\f$. */
  double _thetaMin;
  /** Integration stepsize over \f$\theta\f$. */
  double _ddTheta;
  /** The density of \f$\theta\f$ evaluated on a regular grid from @c _thetaMin with step size @c _ddTheta. */
  Eigen::VectorXd _dtheta;

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


/** The density object for a compound-weibull distribution.
 * This class determines the distribution of the random variable numerically by integrating over
 * the densities of the specified parameter random variables.
 */
class WeibullCompoundDensityObj: public DensityObj
{
protected:
  /** Constructs the density from the given parameter distributions.
   * @param k Specifies the shape parameter random variable.
   * @param lambda Specifies the scale parameter random variable.
   * @param shift Specifies the optional shift, default @c shift=0. */
  WeibullCompoundDensityObj(DensityObj *k, DensityObj *lambda, double shift=0);

public:
  /** Constructs the density from the given parameter distributions.
   * @param k Specifies the shape parameter random variable.
   * @param lambda Specifies the scale parameter random variable.
   * @param shift Specifies the optional shift, default @c shift=0. */
  WeibullCompoundDensityObj(const Var &k, const Var &lambda, double shift=0) throw (AssumptionError);
  /** Destructor. */
  virtual ~WeibullCompoundDensityObj();

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;
  virtual void rangeEst(double alpha, double &a, double &b) const;

protected:
  /** Prepares the numerical itergration on construction. If one of the parameter densities is a
   * delta density, the integration over this parameter is performed "analytically". */
  void _init_int();

protected:
  /** The distribution of the shape parameter. */
  DensityObj *_k;
  /** The distribution of the scale parameter. */
  DensityObj *_lambda;
  /** Shift of the affine transform. */
  double _shift;
  /** Lower bound of the integration interval over \f$k\f$. */
  double _kMin;
  /** Integration stepsize over \f$k\f$. */
  double _ddK;
  /** The density of \f$k\f$ evaluated on a regular grid from @c _kMin with step size @c _ddK. */
  Eigen::VectorXd _dk;
  /** Lower bound of the integration interval over \f$\lambda\f$. */
  double _lambdaMin;
  /** Integration stepsize over \f$\lambda\f$. */
  double _ddLambda;
  /** The density of \f$\lambda\f$ evaluated on a regular grid from @c _lambdaMin with step size
   * @c _ddLambda. */
  Eigen::VectorXd _dlambda;

};


/** Implements a compound-Weibull random variable.
 * That is a Weibull-distributed random variable where the shape and scale parameters are random
 * variables too. */
class WeibullCompoundObj: public CompoundObj
{
public:
  /** Constructs the compound-Weibull random variable from the given parameter random variables. */
  WeibullCompoundObj(const Var &k, const Var &lambda, const std::string &name="");

  void mark();
  virtual Density density();

protected:
  /** A reference to the distribution. */
  WeibullCompoundDensityObj *_density;
};

}

#endif // COMPOUND_HH

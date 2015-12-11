#ifndef __SBB_CONDITIONAL_HH__
#define __SBB_CONDITIONAL_HH__

#include "api.hh"
#include "randomvariable.hh"
#include "density.hh"

namespace stochbb {

/** Implements a conditinal mixture density.
 * That is the density of the random variable
 * \f[
 *  Z = \begin{cases}
 *   Y_1 & \text{ if } X_1 < X_2\\
 *   Y_2 & \text{else.}
 *  \end{cases}
 * \f]
 * assuming \f$X_1,X_2,Y_1\f$ and \fX_1,$X_2,Y_2\f$ mutually independent, the probability of
 * \f$X_1<X_2\f$ is  \f$p_{X_1<X_2}=\intF_{X_1}(x)\,f_{X_2}(x)\,dx\f$ and the density of the
 * mixture is then simply \f$f(x) = p_{X_1<X_2}\,f_{Y_1}(x) + (1- p_{X_1<X_2})f_{Y_2}(x) \f$.
 */
class ConditionalDensityObj: public DensityObj
{
protected:
  ConditionalDensityObj(DensityObj *X1, DensityObj *X2, DensityObj *Y1, DensityObj *Y2);

public:
  ConditionalDensityObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2);

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;

protected:
  DensityObj *_X1;
  DensityObj *_X2;
  DensityObj *_Y1;
  DensityObj *_Y2;
};


/** Implements a conditional mixture.
 * That is
 * \f[
 *  Z = \begin{cases}
 *   Y_1 & \text{ if } X_1 < X_2\\
 *   Y_2 & \text{else.}
 *  \end{cases}
 * \f]
 */
class ConditionalObj: public DerivedVarObj
{
public:
  ConditionalObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2, const std::string &name="");

  virtual void mark();

  virtual Density density();

protected:
  ConditionalDensityObj *_density;
};

}

#endif // __SBB_CONDITIONAL_HH__

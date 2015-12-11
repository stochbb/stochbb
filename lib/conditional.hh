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
 *   Y_1 & \mbox{if }\, X_1 < X_2\\
 *   Y_2 & \mbox{else.}
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


/** Implements the density of a conditional chained variable.
 * That is
 * \f[
 *  Z = \begin{cases}
 *   X_1+Y_1 & \text{ if } X_1 < X_2\\
 *   X_2+Y_2 & \text{else,}
 *  \end{cases}
 * \f]
 * where \f$X_1, X_2, Y_1\f$ and \f$X_1, X_2, Y_2\f$ are mutually independent.
 * This means that Z is the sum (chain) of \f$X_1\f$ and \f$Y_1\f$ if \f$X_1<X_2\f$ and
 * the sum (chain) of \f$X_2\f$ and \f$Y_2\f$ else. Please note that this
 * random variable cannot be implemented useing the @c ConditionalObj class as e.g.
 * \f$X_1+Y_1\f$ (one possible outcome) depends trivially on \f$X_1\f$ (part of the condition).
 */
class CondChainDensityObj: public DensityObj
{
protected:
  CondChainDensityObj(DensityObj *X1, DensityObj *X2, DensityObj *Y1, DensityObj *Y2);

public:
  CondChainDensityObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2);
  virtual void mark();

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  Density affine(double scale, double shift) const;

protected:
  DensityObj *_X1;
  DensityObj *_X2;
  DensityObj *_Y1;
  DensityObj *_Y2;
};

class CondChainObj: public DerivedVarObj
{
public:
  CondChainObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2, const std::string &name="");
  virtual void mark();

  virtual Density density();

protected:
  CondChainDensityObj *_density;
};
}

#endif // __SBB_CONDITIONAL_HH__

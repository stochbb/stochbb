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
 * assuming \f$X_1,X_2,Y_1\f$ and \f$X_1,X_2,Y_2\f$ mutually independent, the probability of
 * \f$X_1<X_2\f$ is  \f$p_{X_1<X_2}=\int F_{X_1}(x)\,f_{X_2}(x)\,dx\f$ and the density of the
 * mixture is then simply \f$f(x) = p_{X_1<X_2}\,f_{Y_1}(x) + (1- p_{X_1<X_2})f_{Y_2}(x) \f$.
 * @ingroup density */
class ConditionalDensityObj: public DensityObj
{
protected:
  /** Constructs a conditional mixture density object from the given density objects. */
  ConditionalDensityObj(DensityObj *X1, DensityObj *X2, DensityObj *Y1, DensityObj *Y2);

public:
  /** Constructs a conditional mixture density object from the given random variables. */
  ConditionalDensityObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2);

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;
  virtual void rangeEst(double alpha, double &a, double &b) const;
  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

protected:
  /** First condition RV. */
  DensityObj *_X1;
  /** Second condition RV. */
  DensityObj *_X2;
  /** First result RV if \f$X_1<X_2\f$. */
  DensityObj *_Y1;
  /** Second result RV if \f$X_1>X_2\f$. */
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
 * @ingroup rv */
class ConditionalObj: public DerivedVarObj
{
public:
  /** Constructs the conditional random variable object from the given random variables. */
  ConditionalObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2,
                 const std::string &name="") throw (Error);

  virtual void mark();

  virtual Density density();

  virtual void sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                      Eigen::Ref<Eigen::MatrixXd> samples) const;

  virtual void print(std::ostream &stream) const;

protected:
  /** A reference to the density object. */
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
 * @ingroup density */
class CondChainDensityObj: public DensityObj
{
protected:
  /** Constructs a @c CondChainDensityObj from the given densities. */
  CondChainDensityObj(DensityObj *X1, DensityObj *X2, DensityObj *Y1, DensityObj *Y2) throw (Error);

public:
  /** Constructs a @c CondChainDensityObj from the given variables. */
  CondChainDensityObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2) throw (Error);
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;
  virtual void rangeEst(double alpha, double &a, double &b) const;

  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &out) const;

protected:
  /** Density of the first condition variable. */
  DensityObj *_X1;
  /** Density of the second condition variable. */
  DensityObj *_X2;
  /** Density of the first result variable. */
  DensityObj *_Y1;
  /** Density of the second result variable. */
  DensityObj *_Y2;
};


/** This class implements a random variable object that represents conditionally chained random
 * variables. That is
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
 * @ingroup rv */
class CondChainObj: public DerivedVarObj
{
public:
  /** Constructs the cond. chained RV object from the given random variables. */
  CondChainObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2, const std::string &name="");
  virtual void mark();

  virtual Density density();

  virtual void print(std::ostream &stream) const;
  virtual void sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                      Eigen::Ref<Eigen::MatrixXd> samples) const;

protected:
  /** Holds a reference to the associated density. */
  CondChainDensityObj *_density;
};

}

#endif // __SBB_CONDITIONAL_HH__

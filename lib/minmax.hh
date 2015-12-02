#ifndef __SBB_MINMAX_HH__
#define __SBB_MINMAX_HH__

#include <vector>
#include "api.hh"
#include "randomvariable.hh"

namespace sbb {

/** Implements the density of a random variable being the maximum of several independent random
 * variables.
 * @ingroup internal */
class MaximumDensityObj: public DensityObj
{
public:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MaximumDensityObj(const std::vector<VarObj *> &variables);
  /** Destructor. */
  virtual ~MaximumDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

  /** Returns a vector of weak references to the underlaying densities. */
  inline const std::vector<DensityObj *> &densities() const { return _densities; }

protected:
  /** The vector of densities. */
  std::vector<DensityObj *> _densities;
};


/** Implements the density of a random variable being the minimum of several independent random
 * variables.
 * @ingroup internal */
class MinimumDensityObj: public DensityObj
{
public:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MinimumDensityObj(const std::vector<VarObj *> &variables);
  /** Destructor. */
  virtual ~MinimumDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

  /** Returns a vector of weak references to the underlaying densities. */
  inline const std::vector<DensityObj *> &densities() const { return _densities; }

protected:
  /** The vector of densities. */
  std::vector<DensityObj *> _densities;
};


/** Implements the random variable being the maximum of several independent random
 * variables.
 * @ingroup internal */
class MaximumObj: public VarObj
{
public:
  /** Constructor from two random variables.
   * @throws AssumptionError If these two random variables are not independent. */
  MaximumObj(const Var &a, const Var &b, const std::string &name="");
  /** Constructor from a vector of random variables.
   * @throws AssumptionError If these random variables are not independent. */
  MaximumObj(const std::vector<Var> &variables, const std::string &name="");
  /** Destructor. */
  virtual ~MaximumObj();
  virtual void mark();

  virtual Density density();

  inline size_t numVariables() const { return _variables.size(); }
  inline Var variable(size_t i) const { _variables[i]->ref(); return _variables[i]; }

protected:
  /** The vector of random variables, this RV depends on. */
  std::vector<VarObj *> _variables;
  /** The density object. */
  MaximumDensityObj *_density;
};


/** Constructs a random variable \f$Y\f$, as the maximum of the given random variables.
 * In contrast to the @c MaximumObj, it handles partial dependency between the random variables
 * if they are formed as sums of independent random variables by first separating the common
 * part. For example
 * \f[
 *  Y = max(X_1+X_2,X_1+X_3) \longrightarrow Y = X_1+max(X_2,X_3)\,,
 * \f]
 * where \f$X_1, X_2\f$ and \f$X_3\f$ are mutually independent random variables.
 * @ingroup internal
 */
VarObj *maximum(const std::vector<VarObj *> &variables);


/** Implements the random variable being the minimum of several independent random
 * variables.
 * @ingroup internal */
class MinimumObj: public VarObj
{
public:
  /** Constructor from two random variables.
   * @throws AssumptionError If these two random variables are not independent. */
  MinimumObj(const Var &a, const Var &b, const std::string &name="");
  /** Constructor from a vector of random variables.
   * @throws AssumptionError If these random variables are not independent. */
  MinimumObj(const std::vector<Var> &variables, const std::string &name="");
  /** Destructor. */
  virtual ~MinimumObj();

  virtual void mark();

  virtual Density density();

  inline size_t numVariables() const { return _variables.size(); }
  inline Var variable(size_t i) const { _variables[i]->ref(); return _variables[i]; }

protected:
  /** The vector of random variables, this RV depends on. */
  std::vector<VarObj *> _variables;
  /** The density object of this random variable. */
  MinimumDensityObj *_density;
};

}

#endif // __SBB_MINMAX_HH__

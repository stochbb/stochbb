#ifndef __SBB_MINMAX_HH__
#define __SBB_MINMAX_HH__

#include <vector>
#include "api.hh"
#include "randomvariable.hh"

namespace sbb {

/** Implements the density of a random variable being the maximum of several independent random
 * variables. */
class MaximumDensityObj: public DensityObj
{
protected:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables (weak references).
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MaximumDensityObj(const std::vector<VarObj *> &variables);

public:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MaximumDensityObj(const std::vector<Var> &variables);

  /** Destructor. */
  virtual ~MaximumDensityObj();

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

  /** Returns a vector of weak references to the underlaying densities. */
  inline const std::vector<DensityObj *> &densities() const { return _densities; }

  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

protected:
  /** The vector of densities. */
  std::vector<DensityObj *> _densities;
  friend class MaximumObj;
};


/** Implements the density of a random variable being the minimum of several independent random
 * variables. */
class MinimumDensityObj: public DensityObj
{
protected:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables (weak references).
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MinimumDensityObj(const std::vector<VarObj *> &variables);

public:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MinimumDensityObj(const std::vector<Var> &variables);

  /** Destructor. */
  virtual ~MinimumDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

  /** Returns a vector of weak references to the underlaying densities. */
  inline const std::vector<DensityObj *> &densities() const { return _densities; }

  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

protected:
  /** The vector of densities. */
  std::vector<DensityObj *> _densities;
  friend class MinimumObj;
};


/** Implements the random variable being the maximum of several independent random
 * variables. */
class MaximumObj: public DerivedVarObj
{
public:
  /** Constructor from a vector of random variables.
   * @throws AssumptionError If these random variables are not independent. */
  MaximumObj(const std::vector<Var> &variables, const std::string &name="");
  /** Destructor. */
  virtual ~MaximumObj();
  virtual void mark();

  virtual Density density();

protected:
  /** The density object. */
  MaximumDensityObj *_density;
};


/** Implements the random variable being the minimum of several independent random
 * variables. */
class MinimumObj: public DerivedVarObj
{
public:
  /** Constructor from a vector of random variables.
   * @throws AssumptionError If these random variables are not independent. */
  MinimumObj(const std::vector<Var> &variables, const std::string &name="");
  /** Destructor. */
  virtual ~MinimumObj();

  virtual void mark();

  virtual Density density();

protected:
  /** The density object of this random variable. */
  MinimumDensityObj *_density;
};

}

#endif // __SBB_MINMAX_HH__

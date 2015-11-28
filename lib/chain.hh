#ifndef CHAIN_HH
#define CHAIN_HH

#include <vector>
#include "randomvariable.hh"

namespace sbb {


class ConvolutionDensityObj: public DensityObj
{
public:
  ConvolutionDensityObj(const std::vector<RandomVariableObj *> &variables);
  virtual ~ConvolutionDensityObj();

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

  inline const std::vector<DensityObj *> &densities() const { return _densities; }

protected:
  std::vector<DensityObj *> _densities;
};


class ChainObj : public RandomVariableObj
{
public:
  ChainObj(RandomVariableObj *a, RandomVariableObj *b);
  ChainObj(const std::vector<RandomVariableObj *> &variables);
  virtual ~ChainObj();
  virtual void mark();

  virtual DensityObj *density();

  inline const std::vector<RandomVariableObj *> &variables() const { return _variables; }

protected:
  std::vector<RandomVariableObj *> _variables;
  ConvolutionDensityObj *_density;
};


/** Implements the random variable
 * \f[
 *  Y = \sum_i X_i\,,
 * \f]
 * where \f$X_i\f$ are independent random variables. */
class Chain: public RandomVariable
{
public:
  typedef ChainObj ObjectType;

public:
  Chain(const RandomVariable &a, const RandomVariable &b);

  Chain &operator =(const Chain &other);

  inline size_t numVariables() const { return _chain->variables().size(); }
  inline RandomVariable variable(size_t idx) const { return _chain->variables()[idx]; }

protected:
  ChainObj *_chain;
};

}

#endif // CHAIN_HH

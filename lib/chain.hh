#ifndef CHAIN_HH
#define CHAIN_HH

#include <vector>
#include "randomvariable.hh"

namespace sbb {


class ConvolutionDensityObj: public DensityObj
{
public:
  ConvolutionDensityObj(const std::vector<VarObj *> &variables);
  virtual ~ConvolutionDensityObj();

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

  inline const std::vector<DensityObj *> &densities() const { return _densities; }

protected:
  std::vector<DensityObj *> _densities;
};


class ChainObj : public VarObj
{
public:
  ChainObj(VarObj *a, VarObj *b);
  ChainObj(const std::vector<VarObj *> &variables);
  virtual ~ChainObj();
  virtual void mark();

  virtual DensityObj *density();

  inline const std::vector<VarObj *> &variables() const { return _variables; }

protected:
  std::vector<VarObj *> _variables;
  ConvolutionDensityObj *_density;
};


/** Implements the random variable
 * \f[
 *  Y = \sum_i X_i\,,
 * \f]
 * where \f$X_i\f$ are independent random variables. */
class Chain: public Var
{
public:
  typedef ChainObj ObjectType;

public:
  Chain(const Var &a, const Var &b);

  Chain &operator =(const Chain &other);

  inline size_t numVariables() const { return _chain->variables().size(); }
  inline Var variable(size_t idx) const { return _chain->variables()[idx]; }

protected:
  ChainObj *_chain;
};

}

#endif // CHAIN_HH

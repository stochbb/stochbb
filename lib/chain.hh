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

}

#endif // CHAIN_HH

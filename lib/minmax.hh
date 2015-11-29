#ifndef __SBB_MINMAX_HH__
#define __SBB_MINMAX_HH__

#include <vector>
#include "randomvariable.hh"

namespace sbb {

class MaximumDensityObj: public DensityObj
{
public:
  MaximumDensityObj(const std::vector<VarObj *> &variables);
  virtual ~MaximumDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

  inline const std::vector<DensityObj *> &densities() const { return _densities; }

protected:
  std::vector<DensityObj *> _densities;
};


class MaximumObj: public VarObj
{
public:
  MaximumObj(VarObj *a, VarObj *b);
  MaximumObj(const std::vector<VarObj *> &variables);
  virtual ~MaximumObj();

  virtual void mark();

  virtual DensityObj *density();

  const std::vector<VarObj *> variables() const {
    return _variables;
  }

protected:
  std::vector<VarObj *> _variables;
  MaximumDensityObj *_density;
};


class Maximum: public Var
{
public:
  typedef MaximumObj ObjectType;

public:
  Maximum(MaximumObj *obj);
  Maximum(const Var &a, const Var &b);
  Maximum(const Maximum &other);
  Maximum &operator=(const Maximum &other);

  inline size_t numVariables() const { return _maximum->variables().size(); }
  inline Var variable(size_t idx) { return _maximum->variables()[idx]; }

protected:
  MaximumObj *_maximum;
};

}

#endif // __SBB_MINMAX_HH__

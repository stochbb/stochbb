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

class MinimumDensityObj: public DensityObj
{
public:
  MinimumDensityObj(const std::vector<VarObj *> &variables);
  virtual ~MinimumDensityObj();
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


class MinimumObj: public VarObj
{
public:
  MinimumObj(VarObj *a, VarObj *b);
  MinimumObj(const std::vector<VarObj *> &variables);
  virtual ~MinimumObj();

  virtual void mark();

  virtual DensityObj *density();

  const std::vector<VarObj *> variables() const {
    return _variables;
  }

protected:
  std::vector<VarObj *> _variables;
  MinimumDensityObj *_density;
};

}

#endif // __SBB_MINMAX_HH__

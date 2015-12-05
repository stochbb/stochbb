#ifndef __SBB_MIXTURE_HH__
#define __SBB_MIXTURE_HH__

#include "api.hh"
#include "randomvariable.hh"
#include "density.hh"


namespace sbb {

class MixtureDensityObj: public DensityObj
{
public:
  MixtureDensityObj(const std::vector<double> &weights, const std::vector<VarObj *> &variables);

  virtual void mark();
  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

protected:
  std::vector<double> _weights;
  std::vector<DensityObj *> _densities;
};


class MixtureObj : public DerivedVarObj
{
public:
  MixtureObj(const std::vector<double> &weights, const std::vector<Var> &variables,
             const std::string &name="");

  virtual void mark();
  virtual Density density();

  inline double weight(size_t i) const { return _weights[i]; }

protected:
  std::vector<double> _weights;
  MixtureDensityObj *_density;
};

}

#endif // __SBB_MIXTURE_HH__

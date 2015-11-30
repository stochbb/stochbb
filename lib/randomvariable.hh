#ifndef __SBB_RANDOMVARIABLE_HH__
#define __SBB_RANDOMVARIABLE_HH__

#include "density.hh"
#include <set>


namespace sbb {


class VarObj: public Object
{
protected:
  VarObj();

public:
  virtual ~VarObj();

  virtual void mark();

  virtual DensityObj *density() = 0;

  inline const std::set<VarObj *> &dependencies() const {
    return _dependencies;
  }

  inline bool dependsOn(VarObj *var) const {
    return 0 != _dependencies.count(var);
  }

  inline bool mutuallyIndep(VarObj *var) const {
    if (dependsOn(var)) { return false; }
    std::set<VarObj *>::const_iterator item = var->dependencies().begin();
    for (; item != var->dependencies().end(); item++) {
      if (dependsOn(*item)) { return false; }
    }
    return true;
  }


protected:
  std::set<VarObj *> _dependencies;
};


class GenericVarObj: public VarObj
{
public:
  GenericVarObj(DensityObj *density);
  virtual ~GenericVarObj();
  virtual void mark();

  virtual DensityObj *density();

public:
  static GenericVarObj *delta(double delay);
  static GenericVarObj *unif(double a, double b);
  static GenericVarObj *norm(double mu, double sigma);
  static GenericVarObj *gamma(double k, double theta);

protected:
  DensityObj *_density;
};


}

#endif // __SBB_RANDOMVARIABLE_HH__

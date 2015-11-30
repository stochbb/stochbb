#ifndef __SBB_RANDOMVARIABLE_HH__
#define __SBB_RANDOMVARIABLE_HH__

#include "density.hh"
#include <set>


namespace sbb {


class VarObj: public Object
{
protected:
  VarObj(const std::string &name="");

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

  inline const std::string &name() const { return _name; }
  inline void setName(const std::string &name) { _name = name; }

protected:
  std::set<VarObj *> _dependencies;
  std::string _name;

protected:
  static std::set<std::string> _var_names;
};


class GenericVarObj: public VarObj
{
public:
  GenericVarObj(DensityObj *density, const std::string &name="");
  virtual ~GenericVarObj();
  virtual void mark();

  virtual DensityObj *density();

public:
  static GenericVarObj *delta(double delay, const std::string &name="");
  static GenericVarObj *unif(double a, double b, const std::string &name="");
  static GenericVarObj *norm(double mu, double sigma, const std::string &name="");
  static GenericVarObj *gamma(double k, double theta, const std::string &name="");

protected:
  DensityObj *_density;
};


}

#endif // __SBB_RANDOMVARIABLE_HH__

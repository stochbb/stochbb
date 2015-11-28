#ifndef __SBB_RANDOMVARIABLE_HH__
#define __SBB_RANDOMVARIABLE_HH__

#include "density.hh"
#include <set>


namespace sbb {


class RandomVariableObj: public Object
{
protected:
  RandomVariableObj();

public:
  virtual ~RandomVariableObj();

  virtual void mark();

  virtual DensityObj *density() = 0;

  inline const std::set<RandomVariableObj *> &dependencies() const {
    return _dependencies;
  }

  inline bool dependsOn(RandomVariableObj *var) const {
    return 0 != _dependencies.count(var);
  }

  inline bool mutuallyIndep(RandomVariableObj *var) const {
    if (dependsOn(var)) { return false; }
    std::set<RandomVariableObj *>::const_iterator item = var->dependencies().begin();
    for (; item != var->dependencies().end(); item++) {
      if (dependsOn(*item)) { return false; }
    }
    return true;
  }


protected:
  std::set<RandomVariableObj *> _dependencies;
};


class GenericRandomVariableObj: public RandomVariableObj
{
public:
  GenericRandomVariableObj(DensityObj *density);
  virtual ~GenericRandomVariableObj();
  virtual void mark();

  virtual DensityObj *density();

public:
  static GenericRandomVariableObj *delta(double delay);
  static GenericRandomVariableObj *unif(double a, double b);
  static GenericRandomVariableObj *norm(double mu, double sigma);

protected:
  DensityObj *_density;
};


class RandomVariable: public Container
{
public:
  typedef RandomVariableObj ObjectType;

public:
  RandomVariable(RandomVariableObj *obj);
  RandomVariable(const RandomVariable &other);

  RandomVariable &operator =(const RandomVariable &other);
  inline RandomVariableObj *operator *() const { return _randomVariable; }

  inline Density density() const {
    return _randomVariable->density();
  }

protected:
  RandomVariableObj *_randomVariable;
};


class GenericRandomVariable: public RandomVariable
{
public:
  typedef GenericRandomVariableObj ObjectType;

public:
  GenericRandomVariable(GenericRandomVariableObj *obj);
  GenericRandomVariable(const Density &density);
  GenericRandomVariable(const GenericRandomVariable &other);

  GenericRandomVariable &operator=(const GenericRandomVariable &other);

public:
  inline static GenericRandomVariable delta(double delay) {
    return GenericRandomVariableObj::delta(delay);
  }

  inline static GenericRandomVariable unif(double a, double b) {
    return GenericRandomVariableObj::unif(a,b);
  }

  inline static GenericRandomVariable norm(double mu, double sigma) {
    return GenericRandomVariableObj::norm(mu, sigma);
  }

protected:
  GenericRandomVariableObj *_genericRV;
};

}

#endif // __SBB_RANDOMVARIABLE_HH__

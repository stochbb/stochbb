#ifndef __SBB_RANDOMVARIABLE_HH__
#define __SBB_RANDOMVARIABLE_HH__

#include "density.hh"
#include <set>
#include <vector>


namespace sbb {

/** Base class of all random variable objects.
 * @ingroup internal */
class VarObj: public Object
{
protected:
  /** Hidden constructor.
   * @param name Specifies the optional name of the variable. */
  VarObj(const std::string &name="");

public:
  /** Destructor. */
  virtual ~VarObj();

  virtual void mark();

  /** Retruns the @c DensityObj of the random variable. */
  virtual DensityObj *density() = 0;

  /** Returns the set of random variables, this RV depends on. */
  inline const std::set<VarObj *> &dependencies() const {
    return _dependencies;
  }

  /** Returns @c true if this random variable depends on the specified one. */
  inline bool dependsOn(VarObj *var) const {
    return 0 != _dependencies.count(var);
  }

  /** Returns @c true if this random variable is mutually independent from the given on. */
  inline bool mutuallyIndep(VarObj *var) const {
    if (dependsOn(var)) { return false; }
    std::set<VarObj *>::const_iterator item = var->dependencies().begin();
    for (; item != var->dependencies().end(); item++) {
      if (dependsOn(*item)) { return false; }
    }
    return true;
  }

  /** Returns the optional name of the random variable. */
  inline const std::string &name() const { return _name; }
  /** Sets the name of the random variable. */
  inline void setName(const std::string &name) { _name = name; }

protected:
  /** The random variables, this RV depends on. */
  std::set<VarObj *> _dependencies;
  /** The optional name of the random variable. */
  std::string _name;

protected:
  /** All names for RVs. */
  static std::set<std::string> _var_names;
};


/** Implements a generic RV defined through its density.
 * @ingroup internal */
class GenericVarObj: public VarObj
{
public:
  /** Constructor from density. */
  GenericVarObj(DensityObj *density, const std::string &name="");
  /** Destructor. */
  virtual ~GenericVarObj();

  virtual void mark();
  virtual DensityObj *density();

public:
  /** Constructs a delta distributed random variable. */
  static GenericVarObj *delta(double delay, const std::string &name="");
  /** Constructs a unifor distributed random variable. */
  static GenericVarObj *unif(double a, double b, const std::string &name="");
  /** Constructs a normal distributed random variable. */
  static GenericVarObj *norm(double mu, double sigma, const std::string &name="");
  /** Constructs a gamma distributed random variable. */
  static GenericVarObj *gamma(double k, double theta, const std::string &name="");

protected:
  /** The density object. */
  DensityObj *_density;
};


/** A set of random variables implementing a proper memory management.
 * Additionally, this class provides methods to compute the union and intersect of
 * two sets.
 * @ingroup internal */
class VarSetObj : public Object
{
public:
  /** Iterator type. */
  typedef std::set<VarObj *>::const_iterator iterator;

public:
  /** Empty constructor. */
  VarSetObj();
  /** Constructor from a std::set of random variables. */
  VarSetObj(const std::set<VarObj *> &variables);
  /** Constructor from a vector of random variables. */
  VarSetObj(const std::vector<VarObj *> &variables);
  /** Copy constructor. */
  VarSetObj(const VarSetObj &other);

  virtual void mark();

  /** Returns @c true if the set is empty. */
  inline bool isEmpty() const {
    return 0 == _vars.size();
  }

  /** Returns the size of the set. */
  inline size_t size() const {
    return _vars.size();
  }

  /** Returns @c true if the set contains the given variable. */
  inline bool contains(VarObj *var) const {
    return (0 != _vars.count(var));
  }

  /** Adds a variable to the set. */
  inline void add(VarObj *var) {
    _vars.insert(var);
  }
  /** Removes the given variable from the set. */
  inline void remove(VarObj *var) {
    _vars.erase(var);
  }

  /** Computes the union with the given set. */
  inline VarSetObj *unite(VarSetObj *other) const {
    VarSetObj *u = new VarSetObj(*this);
    VarSetObj::iterator item = other->begin();
    for (; item != other->end(); item++) {
      u->add(*item);
    }
    return u;
  }

  /** Computes the intersection with the given set. */
  inline VarSetObj *intersect(VarSetObj *other) const {
    VarSetObj *sec = new VarSetObj();
    VarSetObj::iterator item = other->begin();
    for (; item != other->end(); item++) {
      if (this->contains(*item)) {
        sec->add(*item);
      }
    }
    return sec;
  }

  /** Computes the difference between this and the given set. */
  inline VarSetObj *difference(VarSetObj *other) const {
    VarSetObj *diff = new VarSetObj(*this);
    iterator item = other->begin();
    for (; item != other->end(); item++) {
      if (diff->contains(*item)) {
        diff->remove(*item);
      }
    }
    return diff;
  }

  /** Returns the iterator pointing at the first element in the set. */
  inline iterator begin() const {
    return _vars.begin();
  }

  /** Returns the iterator pointing right after the last element in the set. */
  inline iterator end() const {
    return _vars.end();
  }

protected:
  /** The set of variables. */
  std::set<VarObj *> _vars;
};

}

#endif // __SBB_RANDOMVARIABLE_HH__

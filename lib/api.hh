/** @defgroup api Programming interface. */

#ifndef __SBB_API_HH__
#define __SBB_API_HH__

#include "exception.hh"
#include "object.hh"
#include "randomvariable.hh"
#include "density.hh"
#include "chain.hh"
#include "minmax.hh"
#include "simulation.hh"

namespace sbb {

/** Base class of all container classes.
 * @ingroup api */
class Container
{
public:
  /** The object type of the container. */
  typedef Object ObjectType;

protected:
  /** Hidden constructor. */
  Container();
  /** Packs the given objects. */
  explicit Container(Object *obj);
  /** Copy constructor. */
  Container(const Container &other);

public:
  /** Destructor. */
  virtual ~Container();

  /** Assignment operator. */
  const Container &operator=(const Container &other);

  /** Returns @c true if the container is empty. */
  bool isNull() const;

  /** Returns @c true if the object can be casted to the given container. */
  template <class T>
  bool is() const {
    return 0 != dynamic_cast<typename T::ObjectType *>(this->_object);
  }

  /** Casts the object to the given container. */
  template <class T>
  T as() const {
    return T(dynamic_cast<typename T::ObjectType *>(this->_object));
  }

protected:
  /** Boxes the given object. */
  void box(Object *obj);

protected:
  /** Holds the object. */
  Object *_object;
};


class Density: public Container
{
public:
  typedef DensityObj ObjectType;

public:
  Density(DensityObj *obj);
  Density(const Density &other);

  Density &operator=(const Density &other);

  inline void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
    _density->eval(Tmin, Tmax, out);
  }
  inline void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
    _density->evalCDF(Tmin, Tmax, out);
  }
  inline void sample(Eigen::VectorXd &out) const {
    _density->sample(out);
  }

  inline DensityObj *operator *() const { return _density; }

protected:
  DensityObj *_density;
};


class Var: public Container
{
public:
  typedef VarObj ObjectType;

public:
  Var(VarObj *obj);
  Var(const Var &other);

  Var &operator =(const Var &other);

  inline VarObj *operator *() const {
    return _randomVariable;
  }

  inline Density density() const {
    return _randomVariable->density();
  }

  inline const std::string &name() const {
    return _randomVariable->name();
  }
  inline void setName(const std::string &name) {
    _randomVariable->setName(name);
  }

protected:
  VarObj *_randomVariable;
};


class GenericVar: public Var
{
public:
  typedef GenericVarObj ObjectType;

public:
  GenericVar(GenericVarObj *obj);
  GenericVar(const Density &density, const std::string &name="");
  GenericVar(const GenericVar &other);

  GenericVar &operator=(const GenericVar &other);

public:
  inline static GenericVar delta(double delay) {
    return GenericVarObj::delta(delay);
  }

  inline static GenericVar unif(double a, double b) {
    return GenericVarObj::unif(a,b);
  }

  inline static GenericVar norm(double mu, double sigma) {
    return GenericVarObj::norm(mu, sigma);
  }

  inline static GenericVar gamma(double k, double theta) {
    return GenericVarObj::gamma(k, theta);
  }

protected:
  GenericVarObj *_genericRV;
};


/** Implements the random variable
 * \f[
 *  Y = \sum_i X_i\,,
 * \f]
 * where \f$X_i\f$ are independent random variables.
 * @ingroup api */
class Chain: public Var
{
public:
  typedef ChainObj ObjectType;

public:
  Chain(const Var &a, const Var &b, const std::string &name="");

  Chain &operator =(const Chain &other);

  inline size_t numVariables() const { return _chain->variables().size(); }
  inline Var variable(size_t idx) const { return _chain->variables()[idx]; }

protected:
  ChainObj *_chain;
};


class Maximum: public Var
{
public:
  typedef MaximumObj ObjectType;

public:
  Maximum(MaximumObj *obj);
  Maximum(const Var &a, const Var &b, const std::string &name="");
  Maximum(const Maximum &other);
  Maximum &operator=(const Maximum &other);

  inline size_t numVariables() const { return _maximum->variables().size(); }
  inline Var variable(size_t idx) { return _maximum->variables()[idx]; }

protected:
  MaximumObj *_maximum;
};


class Minimum: public Var
{
public:
  typedef MinimumObj ObjectType;

public:
  Minimum(MinimumObj *obj);
  Minimum(const Var &a, const Var &b);
  Minimum(const Minimum &other);
  Minimum &operator=(const Minimum &other);

  inline size_t numVariables() const { return _minimum->variables().size(); }
  inline Var variable(size_t idx) { return _minimum->variables()[idx]; }

protected:
  MinimumObj *_minimum;
};


class Simulation: public Container
{
public:
  typedef SimulationObj ObjectType;

public:
  Simulation();
  Simulation(SimulationObj *object);
  Simulation(const Simulation &other);

  Simulation &operator= (const Simulation &other);

  inline bool hasVar(const std::string &id) const {
    return _simulation->hasVar(id);
  }
  inline Var var(const std::string &id) const {
    return _simulation->var(id);
  }
  inline void addVar(const std::string &id, Var &var) const {
    _simulation->addVar(id, *var);
  }

  inline double tMin() const { return _simulation->tMin(); }
  inline void setTMin(double tMin) { _simulation->setTMin(tMin); }
  inline double tMax() const { return _simulation->tMax(); }
  inline void setTMax(double tMax) { _simulation->setTMax(tMax); }
  inline size_t steps() const { return _simulation->steps(); }
  inline void setSteps(size_t steps) const { _simulation->setSteps(steps); }

  inline size_t numOutputVars() const { return _simulation->outputVars().size(); }
  inline Var outputVar(size_t idx) const { return _simulation->outputVars()[idx]; }
  inline void addOutputVar(const Var &var) { return _simulation->addOutputVar(*var); }

  inline void run(Eigen::MatrixXd &out) const { return _simulation->run(out); }

public:
  static Simulation fromXml(const std::string &filename);

protected:
  SimulationObj *_simulation;
};


}

#include "operators.hh"
#include "logger.hh"

#endif // __SBB_API_HH__

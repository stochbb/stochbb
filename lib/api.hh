/** @defgroup api User and application programming interface
 * This group collects all classes provided as an API to the core functionality of libstochbb.
 * Instead of using the @ref internal classes, the user should resort to API classes as they
 * ease the memory management a lot.
 *
 * All API classes are derived from the @c Container class which keep track of all objects being
 * directly accessible. The @c GC class, then implements a mark and sweep garbage collector,
 * freeing all objects that are not reachable anymore. Although the container classes are passed
 * around like values, they only consist of a reference to the actual @c Object instance. Hence
 * it is a very efficient way of handling pointers or references without the burden of keeping
 * track of unreacable memory. */

#ifndef __SBB_API_HH__
#define __SBB_API_HH__

#include <vector>
#include <Eigen/Eigen>
#include "exception.hh"
#include "object.hh"


namespace sbb {

// Forward declarations
class DensityObj;
class AtomicDensityObj;
class VarObj;
class AtomicVarObj;
class DerivedVarObj;
class ChainObj;
class MaximumObj;
class MinimumObj;
class SimulationObj;


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
  /** Packs the given objects. Takes the reference of the object. */
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

  /** Returns a weak reference to the object. */
  Object *operator *() const { return _object; }
  /** Returns a weak reference to the object. */
  Object *operator ->() const { return _object; }
  /** Returns a new reference to the object. */
  Object *ref() const { _object->ref(); return _object; }

  /** Returns @c true if the object can be casted to the given container. */
  template <class T>
  bool is() const {
    return 0 != dynamic_cast<typename T::ObjectType *>(this->_object);
  }

  /** Casts the object to the given container. */
  template <class T>
  T as() const {
    typename T::ObjectType *obj = dynamic_cast<typename T::ObjectType *>(this->_object);
    if (obj) { obj->ref(); }
    return T(obj);
  }

protected:
  /** Holds the object. */
  Object *_object;
};


/** Base class of all densities (PDFs).
 * @ingroup api */
class Density: public Container
{
public:
  /** The object type of the container. */
  typedef DensityObj ObjectType;

public:
  /** Packs the given @c DensityObj and taks the reference. */
  Density(DensityObj *obj);
  /** Copy constructor. */
  Density(const Density &other);
  /** Assigment operator. */
  Density &operator=(const Density &other);

  /** Evaluates the density (PDF) on a regular grid \f$[Tmin, Tmax)\f$ where the number
   * of grid points is specified via the length of the output vector @c out. The results are
   * stored into the output vector. */
  void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;

  /** Evaluates the probability function (CDF) on a regular grid \f$[Tmin, Tmax)\f$ where the number
   * of grid points is specified via the length of the output vector @c out. The results are
   * stored into the output vector. */
  void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

  /** Retruns a weak reference to the @c DensityObj. */
  inline DensityObj *operator *() const { return _density; }
  /** Retruns a weak reference to the @c DensityObj. */
  inline DensityObj *operator ->() const { return _density; }
  /** Retruns a new reference to the @c DensityObj. */
  inline DensityObj *ref() const { _object->ref(); return _density; }

protected:
  /** Holds the density object. */
  DensityObj *_density;
};


/** The densities of @c AtomicVar.
 * This is the base class of all densities of atomic random variables, that is a random variable
 * that does not depend on others.
 * @ingroup api */
class AtomicDensity: public Density
{
public:
  /** The object type of the container. */
  typedef AtomicDensityObj ObjectType;

public:
  /** Packs the density object and takes the reference. */
  AtomicDensity(AtomicDensityObj *obj);
  /** Copy constructor. */
  AtomicDensity(const AtomicDensity &other);
  /** Assignment operator. */
  AtomicDensity &operator =(const AtomicDensity &other);

  /** Retruns a weak reference to the @c AtomicDensityObj. */
  inline AtomicDensityObj *operator *() const { return _atomic_density; }
  /** Retruns a weak reference to the @c AtomicDensityObj. */
  inline AtomicDensityObj *operator ->() const { return _atomic_density; }
  /** Retruns a new reference to the @c AtomicDensityObj. */
  inline AtomicDensityObj *ref() const { _object->ref(); return _atomic_density; }

  /** Samples from the random variable. */
  void sample(Eigen::VectorXd &out) const;

protected:
  /** Holds the reference to the object instance. */
  AtomicDensityObj *_atomic_density;
};


/** Base class of all random variables.
 * @ingroup api */
class Var: public Container
{
public:
  /** The object type of the container. */
  typedef VarObj ObjectType;

public:
  /** Empty constructor. */
  Var();
  /** Packs the given random variable object and takes the reference. */
  Var(VarObj *obj);
  /** Copy constructor. */
  Var(const Var &other);
  /** Assingment operator. */
  Var &operator =(const Var &other);

  /** Returns a weak reference to the random variable object. */
  inline VarObj *operator *() const {
    return _randomVariable;
  }
  /** Returns a weak reference to the random variable object. */
  inline VarObj *operator ->() const {
    return _randomVariable;
  }
  /** Returns a new reference to the random variable object. */
  inline VarObj *ref() const {
    _object->ref(); return _randomVariable;
  }

  /** Retruns true if the given variable container refers to the same random variable. */
  inline bool operator==(const Var &other) const {
        return _randomVariable == other._randomVariable;
  }
  /** Implements a partial ordering of random variables. */
  inline bool operator<(const Var &other) const {
    return _randomVariable < other._randomVariable;
  }

  /** Returns a reference to the density associated with this random variable. */
  Density density() const;

  /** Returns @c true if this random variable depends on the given one. */
  bool dependsOn(const Var &other) const;
  /** Returns @c true if this and the given random variable are mutually independent. */
  bool mutuallyIndep(const Var &other) const;

  /** Returns the optional name of the random variable. */
  const std::string &name() const;
  /** Sets the name of the random variable. */
  void setName(const std::string &name);

protected:
  /** Holds the pointer to the random variable object. */
  VarObj *_randomVariable;
};


/** A generic random variable derived from a chosen density.
 * @ingroup api */
class AtomicVar: public Var
{
public:
  /** The object type of the container. */
  typedef AtomicVarObj ObjectType;

public:
  /** Packs the given @c GenericVarObj. */
  AtomicVar(AtomicVarObj *obj);
  /** Constructs a generic random variable from distribution. */
  AtomicVar(const AtomicDensity &density, const std::string &name="");
  /** Copy constructor. */
  AtomicVar(const AtomicVar &other);
  /** Assignement operator. */
  AtomicVar &operator=(const AtomicVar &other);

public:
  /** Constructs a delta-distributed "random" variable. */
  static AtomicVar delta(double delay);
  /** Constructs a uniformly distributed random variable. */
  static AtomicVar unif(double a, double b);
  /** Constructs a normal distributed random variable. */
  static AtomicVar norm(double mu, double sigma);
  /** Constructs a gamma distributed random variable. */
  static AtomicVar gamma(double k, double theta);

protected:
  /** Holds the @c GenericVarObj. */
  AtomicVarObj *_genericRV;
};


class DerivedVar: public Var
{
public:
  typedef DerivedVarObj ObjectType;

protected:
  DerivedVar(DerivedVarObj *obj);

public:
  DerivedVar(const DerivedVar &other);
  DerivedVar &operator =(const DerivedVar &other);

  /** Returns the number of random variables the variable depends on directly. */
  size_t numVariables() const;
  /** Returns the i-th random variable the chain depends on directly. */
  Var variable(size_t idx) const;

protected:
  DerivedVarObj *_derived_var;
};


/** Implements the sum of the independent random variables \f$X_i, i=1,\dots,N\f$,
 *  \f$Y = \sum_i X_i\f$.
 *
 * This class is called @c Chain as it implements a chain of random processes. Assume the
 * independent random processes \f$X_1\f$ and \f$X_2\f$ and \f$X_1\f$ triggers \f$X_2\f$,
 * than the time of the second process (\f$X_2\f$) to complete is simply the sum of times
 * of both processes.
 *
 * @ingroup api */
class Chain: public DerivedVar
{
public:
  /** The object type of the container. */
  typedef ChainObj ObjectType;

public:
  /** Packs the given @c ChainObj. */
  Chain(ChainObj *obj);
  /** Constructs a chain of two random variables. */
  Chain(const std::vector<Var> &variables, const std::string &name="");
  /** Assignment operator. */
  Chain &operator =(const Chain &other);

protected:
  /** Holds the @c ChainObj. */
  ChainObj *_chain;
};


/** Implements the maximum of the \f$N\f$ independent random variables \f$X_i,i=1,\dots,N\f$,
 * \f$Y = \max\{X_1,\dots, X_N\}\f$.
 *
 * In terms of random processes: Consider \f$N\f$ parallel processes \f$X_i,i=1,\dots,N\f$, then
 * this class will represent the random variable of the waiting time of the last process.
 *
 * @ingroup api */
class Maximum: public DerivedVar
{
public:
  /** Object type of the container. */
  typedef MaximumObj ObjectType;

public:
  /** Packs the given @c MaximumObj. */
  Maximum(MaximumObj *obj);
  /** Constructs a maximum random variable from the given ones. */
  Maximum(const std::vector<Var> &variables, const std::string &name="");
  /** Copy constructor. */
  Maximum(const Maximum &other);
  /** Assignment operator. */
  Maximum &operator=(const Maximum &other);

protected:
  /** Holds the @c MaximumObj. */
  MaximumObj *_maximum;
};


/** Implements the minimum of the \f$N\f$ independent random variables \f$X_i\f$,
 * \f$Y = \min\{X_1,\dots, X_n\}\f$.
 *
 * In terms of random processes: Consider \f$N\f$ parallel processes \f$X_i,i=1,\dots,N\f$, then
 * this class will represent the random variable of the waiting time of the first process completed.
 *
 * @ingroup api */
class Minimum: public DerivedVar
{
public:
  /** The object type of the container. */
  typedef MinimumObj ObjectType;

public:
  /** Packs the given @c MinimumObj. */
  Minimum(MinimumObj *obj);
  /** Constructs a minimum random variable from the given ones. */
  Minimum(const std::vector<Var> &variables);
  /** Copy constructor. */
  Minimum(const Minimum &other);
  /** Assignment operator. */
  Minimum &operator=(const Minimum &other);

protected:
  /** Holds the @c MinimumObj. */
  MinimumObj *_minimum;
};


/** Collects several variable definitions and which PDFs are evaluated.
 *
 * The simplest way to construct a @c Simulation is to parse a simulation specification from xml
 * e.g.
 * @code
 *  Simulation sim = Simulation::fromXml("filename.xml");
 * @endcode
 *
 * @ingroup api */
class Simulation: public Container
{
public:
  /** The object type of the container. */
  typedef SimulationObj ObjectType;

public:
  /** Empty constructor. */
  Simulation();
  /** Packs the given @c SimulationObj. */
  Simulation(SimulationObj *object);
  /** Copy constructor. */
  Simulation(const Simulation &other);
  /** Assignment operator. */
  Simulation &operator= (const Simulation &other);

  /** Returns @c true if the given random variable is known. */
  bool hasVar(const std::string &id) const;
  /** Returns the specified random variable. */
  Var var(const std::string &id) const;
  /** Adds a random variable to the simulation. */
  void addVar(const std::string &id, Var &var) const;

  /** Returns the start time of the simulation. */
  double tMin() const;
  /** Sets the start time of the simulation. */
  void setTMin(double tMin);
  /** Returns the end time of the simulation. */
  double tMax() const;
  /** Sets the end time of the simulation. */
  void setTMax(double tMax);
  /** Returns the number of time-steps. */
  size_t steps() const;
  /** Sets the number of time-steps. */
  void setSteps(size_t steps) const;
  /** Returns the number of output variables. */
  size_t numOutputVars() const;
  /** Returns the specified output variable. */
  Var outputVar(size_t idx) const;
  /** Adds a output variable to the simulation. */
  void addOutputVar(const Var &var);

  /** Performs the simulation and stores the results into the given matrix.
   * The matrix gets resized. Each column represents a output variable where the first
   * column is time. Each row represetns a time-point at which the PDFs are evaluated. */
  void run(Eigen::MatrixXd &out) const;

public:
  /** Parses the simulation specification from XML.
   * @returns A Simulation instance.
   * @throws ParserError If parsing fails. */
  static Simulation fromXml(const std::string &filename);

protected:
  /** Holds the SimulationObj instance. */
  SimulationObj *_simulation;
};


}

#include "operators.hh"
#include "logger.hh"

#endif // __SBB_API_HH__

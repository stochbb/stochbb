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


/** Base class of all densities (PDFs).
 * @ingroup api */
class Density: public Container
{
public:
  /** The object type of the container. */
  typedef DensityObj ObjectType;

public:
  /** Packs the given @c DensityObj. */
  Density(DensityObj *obj);
  /** Copy constructor. */
  Density(const Density &other);
  /** Assigment operator. */
  Density &operator=(const Density &other);

  /** Evaluates the density (PDF) on a regular grid \f$[Tmin, Tmax)\f$ where the number
   * of grid points is specified via the length of the output vector @c out. The results are
   * stored into the output vector. */
  inline void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
    _density->eval(Tmin, Tmax, out);
  }
  /** Evaluates the probability function (CDF) on a regular grid \f$[Tmin, Tmax)\f$ where the number
   * of grid points is specified via the length of the output vector @c out. The results are
   * stored into the output vector. */
  inline void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
    _density->evalCDF(Tmin, Tmax, out);
  }
  /** Samples from the distribution. */
  inline void sample(Eigen::VectorXd &out) const {
    _density->sample(out);
  }
  /** Retruns a reference to the @c DensityObj. */
  inline DensityObj *operator *() const { return _density; }

protected:
  /** Holds the density object. */
  DensityObj *_density;
};


/** Base class of all random variables.
 * @ingroup api */
class Var: public Container
{
public:
  /** The object type of the container. */
  typedef VarObj ObjectType;

public:
  /** Packs the given random variable object. */
  Var(VarObj *obj);
  /** Copy constructor. */
  Var(const Var &other);
  /** Assingment operator. */
  Var &operator =(const Var &other);
  /** Returns a pointer to the random variable object. */
  inline VarObj *operator *() const {
    return _randomVariable;
  }

  /** Returns a reference to the density associated with this random variable. */
  inline Density density() const {
    return _randomVariable->density();
  }
  /** Returns the optional name of the random variable. */
  inline const std::string &name() const {
    return _randomVariable->name();
  }
  /** Sets the name of the random variable. */
  inline void setName(const std::string &name) {
    _randomVariable->setName(name);
  }

protected:
  /** Holds the pointer to the random variable object. */
  VarObj *_randomVariable;
};


/** A generic random variable derived from a chosen density.
 * @ingroup api */
class GenericVar: public Var
{
public:
  /** The object type of the container. */
  typedef GenericVarObj ObjectType;

public:
  /** Packs the given @c GenericVarObj. */
  GenericVar(GenericVarObj *obj);
  /** Constructs a generic random variable from distribution. */
  GenericVar(const Density &density, const std::string &name="");
  /** Copy constructor. */
  GenericVar(const GenericVar &other);
  /** Assignement operator. */
  GenericVar &operator=(const GenericVar &other);

public:
  /** Constructs a delta-distributed "random" variable. */
  inline static GenericVar delta(double delay) {
    return GenericVarObj::delta(delay);
  }
  /** Constructs a uniformly distributed random variable. */
  inline static GenericVar unif(double a, double b) {
    return GenericVarObj::unif(a,b);
  }
  /** Constructs a normal distributed random variable. */
  inline static GenericVar norm(double mu, double sigma) {
    return GenericVarObj::norm(mu, sigma);
  }
  /** Constructs a gamma distributed random variable. */
  inline static GenericVar gamma(double k, double theta) {
    return GenericVarObj::gamma(k, theta);
  }

protected:
  /** Holds the @c GenericVarObj. */
  GenericVarObj *_genericRV;
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
class Chain: public Var
{
public:
  /** The object type of the container. */
  typedef ChainObj ObjectType;

public:
  /** Packs the given @c ChainObj. */
  Chain(ChainObj *obj);
  /** Constructs a chain of two random variables. */
  Chain(const Var &a, const Var &b, const std::string &name="");
  /** Assignment operator. */
  Chain &operator =(const Chain &other);

  /** Returns the number of random variables the chain depends on. */
  inline size_t numVariables() const { return _chain->variables().size(); }
  /** Returns the i-th random variable the chain depends on. */
  inline Var variable(size_t idx) const { return _chain->variables()[idx]; }

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
class Maximum: public Var
{
public:
  /** Object type of the container. */
  typedef MaximumObj ObjectType;

public:
  /** Packs the given @c MaximumObj. */
  Maximum(MaximumObj *obj);
  /** Constructs a maximum random variable from the given ones. */
  Maximum(const Var &a, const Var &b, const std::string &name="");
  /** Copy constructor. */
  Maximum(const Maximum &other);
  /** Assignment operator. */
  Maximum &operator=(const Maximum &other);

  /** Returns the number of random variables the maximum depends on. */
  inline size_t numVariables() const { return _maximum->variables().size(); }
  /** Returns the i-th random variable the maximum depends on. */
  inline Var variable(size_t idx) { return _maximum->variables()[idx]; }

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
class Minimum: public Var
{
public:
  /** The object type of the container. */
  typedef MinimumObj ObjectType;

public:
  /** Packs the given @c MinimumObj. */
  Minimum(MinimumObj *obj);
  /** Constructs a minimum random variable from the given ones. */
  Minimum(const Var &a, const Var &b);
  /** Copy constructor. */
  Minimum(const Minimum &other);
  /** Assignment operator. */
  Minimum &operator=(const Minimum &other);

  /** Returns the number of random variables the minimum depends on. */
  inline size_t numVariables() const { return _minimum->variables().size(); }
  /** Returns the i-th random variable the minimum depends on. */
  inline Var variable(size_t idx) { return _minimum->variables()[idx]; }

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
  inline bool hasVar(const std::string &id) const {
    return _simulation->hasVar(id);
  }
  /** Returns the specified random variable. */
  inline Var var(const std::string &id) const {
    return _simulation->var(id);
  }
  /** Adds a random variable to the simulation. */
  inline void addVar(const std::string &id, Var &var) const {
    _simulation->addVar(id, *var);
  }

  /** Returns the start time of the simulation. */
  inline double tMin() const { return _simulation->tMin(); }
  /** Sets the start time of the simulation. */
  inline void setTMin(double tMin) { _simulation->setTMin(tMin); }
  /** Returns the end time of the simulation. */
  inline double tMax() const { return _simulation->tMax(); }
  /** Sets the end time of the simulation. */
  inline void setTMax(double tMax) { _simulation->setTMax(tMax); }
  /** Returns the number of time-steps. */
  inline size_t steps() const { return _simulation->steps(); }
  /** Sets the number of time-steps. */
  inline void setSteps(size_t steps) const { _simulation->setSteps(steps); }
  /** Returns the number of output variables. */
  inline size_t numOutputVars() const { return _simulation->outputVars().size(); }
  /** Returns the specified output variable. */
  inline Var outputVar(size_t idx) const { return _simulation->outputVars()[idx]; }
  /** Adds a output variable to the simulation. */
  inline void addOutputVar(const Var &var) { return _simulation->addOutputVar(*var); }

  /** Performs the simulation and stores the results into the given matrix.
   * The matrix gets resized. Each column represents a output variable where the first
   * column is time. Each row represetns a time-point at which the PDFs are evaluated. */
  inline void run(Eigen::MatrixXd &out) const { return _simulation->run(out); }

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

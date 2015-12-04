/** \page api Application programming interface
 * This group collects all classes provided as an API to the core functionality of libstochbb.
 * Instead of using the internal classes, the user should resort to API classes as they
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
class ExactSamplerObj;
class MarginalSamplerObj;


/** Base class of all container classes. */
class Container
{
public:
  /** @cond internal
   * The object type of the container. */
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

  /** Returns a weak reference to the object. */
  Object *operator *() const { return _object; }
  /** Returns a weak reference to the object. */
  Object *operator ->() const { return _object; }
  /** Returns a new reference to the object. */
  Object *ref() const { _object->ref(); return _object; }
  //! @endcond

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
    typename T::ObjectType *obj = dynamic_cast<typename T::ObjectType *>(this->_object);
    if (obj) { obj->ref(); }
    return T(obj);
  }

protected:
  /** @cond internal
   * Holds the object. */
  Object *_object;
  //! @endcond

};


/** Base class of all densities (PDFs). */
class Density: public Container
{
public:
  /** @cond internal
   * The object type of the container. */
  typedef DensityObj ObjectType;

public:
  /** Packs the given @c DensityObj and taks the reference. */
  Density(DensityObj *obj);
  /** Copy constructor. */
  Density(const Density &other);
  /** Assigment operator. */
  Density &operator=(const Density &other);
  //! @endcond

  /** Evaluates the density (PDF) on a regular grid \f$[Tmin, Tmax)\f$ where the number
   * of grid points is specified via the length of the output vector @c out. The results are
   * stored into the output vector. */
  void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;

  /** Evaluates the probability function (CDF) on a regular grid \f$[Tmin, Tmax)\f$ where the number
   * of grid points is specified via the length of the output vector @c out. The results are
   * stored into the output vector. */
  void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

  //! @cond internal
  /** Retruns a weak reference to the @c DensityObj. */
  inline DensityObj *operator *() const { return _density; }
  /** Retruns a weak reference to the @c DensityObj. */
  inline DensityObj *operator ->() const { return _density; }
  /** Retruns a new reference to the @c DensityObj. */
  inline DensityObj *ref() const { _object->ref(); return _density; }

protected:
  /** Holds the density object. */
  DensityObj *_density;
  //! @endcond
};


/** The densities of @c AtomicVar.
 * This is the base class of all densities of atomic random variables, that is a random variable
 * that does not depend on others. */
class AtomicDensity: public Density
{
public:
  //! @cond internal
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
  //! @endcond

  /** Samples from the random variable. */
  void sample(Eigen::VectorXd &out) const;

protected:
  //! @cond internal
  /** Holds the reference to the object instance. */
  AtomicDensityObj *_atomic_density;
  //! @endcond
};


/** Base class of all random variables. */
class Var: public Container
{
public:
  //! @cond internal
  /** The object type of the container. */
  typedef VarObj ObjectType;
  //! @endcond
public:
  /** Empty constructor. */
  Var();

  //! @cond internal
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
  //! @endcond

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
  //! @cond internal
  /** Holds the pointer to the random variable object. */
  VarObj *_randomVariable;
  //! @endcond
};


/** A generic random variable derived from a chosen density. */
class AtomicVar: public Var
{
public:
  //! @cond internal
  /** The object type of the container. */
  typedef AtomicVarObj ObjectType;

public:
  /** Packs the given @c GenericVarObj. */
  AtomicVar(AtomicVarObj *obj);
  //! @endcond

  /** Constructs a generic random variable from distribution. */
  AtomicVar(const AtomicDensity &density, const std::string &name="");

  //! @cond internal
  /** Copy constructor. */
  AtomicVar(const AtomicVar &other);
  /** Assignement operator. */
  AtomicVar &operator=(const AtomicVar &other);
  //! @endcond

  /** Constructs a delta-distributed "random" variable. */
  static AtomicVar delta(double delay);
  /** Constructs a uniformly distributed random variable. */
  static AtomicVar unif(double a, double b);
  /** Constructs a normal distributed random variable. */
  static AtomicVar norm(double mu, double sigma);
  /** Constructs a gamma distributed random variable. */
  static AtomicVar gamma(double k, double theta);

protected:
  //! @cond internal
  /** Holds the @c GenericVarObj. */
  AtomicVarObj *_genericRV;
  //! @endcond
};


/** Base class of all derived random variables. */
class DerivedVar: public Var
{
public:
  //! @cond internal
  /** The object type of the container. */
  typedef DerivedVarObj ObjectType;

protected:
  /** Hidden constructor. */
  DerivedVar(DerivedVarObj *obj);

public:
  /** Copy constructor. */
  DerivedVar(const DerivedVar &other);
  /** Assignment operator. */
  DerivedVar &operator =(const DerivedVar &other);
  //! @endcond

  /** Returns the number of random variables the variable depends on directly. */
  size_t numVariables() const;
  /** Returns the i-th random variable the chain depends on directly. */
  Var variable(size_t idx) const;

protected:
  //! @cond internal
  /** Holds a reference to the @c DerivedVarObj instance. */
  DerivedVarObj *_derived_var;
  //! @endcond
};


/** Implements the sum of the independent random variables \f$X_i, i=1,\dots,N\f$,
 *  \f$Y = \sum_i X_i\f$.
 *
 * This class is called @c Chain as it implements a chain of random processes. Assume the
 * independent random processes \f$X_1\f$ and \f$X_2\f$ and \f$X_1\f$ triggers \f$X_2\f$,
 * than the time of the second process (\f$X_2\f$) to complete is simply the sum of times
 * of both processes. */
class Chain: public DerivedVar
{
public:
  //! @cond internal
  /** The object type of the container. */
  typedef ChainObj ObjectType;

public:
  /** Packs the given @c ChainObj. */
  Chain(ChainObj *obj);
  //! @endcond

  /** Constructs a chain of two random variables. */
  Chain(const std::vector<Var> &variables, const std::string &name="");

  //! @cond internal
  /** Assignment operator. */
  Chain &operator =(const Chain &other);

protected:
  /** Holds the @c ChainObj. */
  ChainObj *_chain;
  //! @endcond
};


/** Implements the maximum of the \f$N\f$ independent random variables \f$X_i,i=1,\dots,N\f$,
 * \f$Y = \max\{X_1,\dots, X_N\}\f$.
 *
 * In terms of random processes: Consider \f$N\f$ parallel processes \f$X_i,i=1,\dots,N\f$, then
 * this class will represent the random variable of the waiting time of the last process. */
class Maximum: public DerivedVar
{
public:
  //! @cond internal
  /** Object type of the container. */
  typedef MaximumObj ObjectType;

public:
  /** Packs the given @c MaximumObj. */
  Maximum(MaximumObj *obj);
  //! @endcond

  /** Constructs a maximum random variable from the given ones. */
  Maximum(const std::vector<Var> &variables, const std::string &name="");

  //! @cond internal
  /** Copy constructor. */
  Maximum(const Maximum &other);
  /** Assignment operator. */
  Maximum &operator=(const Maximum &other);

protected:
  /** Holds the @c MaximumObj. */
  MaximumObj *_maximum;
  //! @endcond
};


/** Implements the minimum of the \f$N\f$ independent random variables \f$X_i\f$,
 * \f$Y = \min\{X_1,\dots, X_n\}\f$.
 *
 * In terms of random processes: Consider \f$N\f$ parallel processes \f$X_i,i=1,\dots,N\f$, then
 * this class will represent the random variable of the waiting time of the first process completed. */
class Minimum: public DerivedVar
{
public:
  //! @cond internal
  /** The object type of the container. */
  typedef MinimumObj ObjectType;

public:
  /** Packs the given @c MinimumObj. */
  Minimum(MinimumObj *obj);
  //! @endcond

  /** Constructs a minimum random variable from the given ones. */
  Minimum(const std::vector<Var> &variables);

  //! @cond internal
  /** Copy constructor. */
  Minimum(const Minimum &other);
  /** Assignment operator. */
  Minimum &operator=(const Minimum &other);

protected:
  /** Holds the @c MinimumObj. */
  MinimumObj *_minimum;
  //! @endcond
};


/** Collects several variable definitions and which PDFs are evaluated.
 *
 * The simplest way to construct a @c Simulation is to parse a simulation specification from xml
 * e.g.
 * @code
 *  Simulation sim = Simulation::fromXml("filename.xml");
 * @endcode
 */
class Simulation: public Container
{
public:
  //! @cond internal
  /** The object type of the container. */
  typedef SimulationObj ObjectType;

public:
  /** Packs the given @c SimulationObj. */
  Simulation(SimulationObj *object);
  /** Copy constructor. */
  Simulation(const Simulation &other);
  /** Assignment operator. */
  Simulation &operator= (const Simulation &other);
  //! @endcond

  /** Empty constructor. */
  Simulation();

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
  //! @cond internal
  /** Holds the SimulationObj instance. */
  SimulationObj *_simulation;
  //! @endcond
};


/** Implements a sampler for several possibly dependent random variables.
 * This class should be used if one is interesed in sampling several random variables
 * simultaneously. If only one variable of a very large system of random variables is sampled,
 * consider also the @c MarginalSampler class. */
class ExactSampler: public Container
{
public:
  //! @cond internal
  /** Object type of the container. */
  typedef ExactSamplerObj ObjectType;
  //! @endcond

public:
  /** Constructs a sampler for the given random variables. */
  ExactSampler(const std::vector<Var> &variables);

  //! @cond internal
  /** Copy constructor. */
  ExactSampler(const ExactSampler &other);
  /** Assignment operator. */
  ExactSampler &operator =(const ExactSampler &other);
  //! @endcond

  /** Samples from the random variables passed to the constructor.
   * The number of row of @c out specify the number of samples to draw. The number of columns
   * must match the number of random variables passed to the constructor, i.e. each colum
   * corresponds to a random variable passed to the constructor. */
  void sample(Eigen::MatrixXd &out) const;

protected:
  //! @cond internal
  /** A reference to the sampler object. */
  ExactSamplerObj *_sampler;
  //! @endcond
};


/** Implements an approximative marginal sampler for a sinlge random variable.
 * This sampler first evaluates the CDF of the random variable passed to the constructor either
 * numerically or if possible analytically on a specified grid. Then a linear interpolation of
 * its inverse is used to generate samples which are distributed approximately like the chosen
 * random variable.
 *
 * This sampler can be used if one wants to sample from a single variable of a large system of
 * dependent random variables where the exact sampling method (see @c ExactSampler) is to slow.
 * This sampler, however, only samples from a marginal. Hence the @c ExactSampler must be used if
 * more than one mutually dependent variables are sampled simultaneously. */
class MarginalSampler: public Container
{
public:
  //! @cond internal
  /** Object type of the container. */
  typedef MarginalSamplerObj ObjectType;
  //! @endcond

public:
  /** Constructs a sampler for the specified variable.
   * @param var Specifies the random variable to sample from.
   * @param Tmin Specifies the lower-bound of the grid on which the CDF is evaluated.
   * @param Tmax Specifies the upper-bound of the grid on which the CDF is evaluated.
   * @param steps Specifies the number of grid points. */
  MarginalSampler(const Var &var, double Tmin, double Tmax, size_t steps);

  //! @cond internal
  /** Copy constructor. */
  MarginalSampler(const MarginalSampler &other);
  /** Assignment operator. */
  MarginalSampler &operator =(const MarginalSampler &other);
  //! @endcond

  /** Samples from the marginal. The number of samples is specified by the number of elements
   * in the vector @c out. */
  void sample(Eigen::VectorXd &out) const;

protected:
  //! @cond internal
  /** Holds a reference to the sampler object. */
  MarginalSamplerObj *_sampler;
  //! @endcond
};

}

#include "operators.hh"
#include "logger.hh"

#endif // __SBB_API_HH__

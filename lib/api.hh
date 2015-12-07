/** \page api Application programming interface
 * The application programming interface (API) allows to assemble processes programmatically in C++.
 * All API classes are derived from the @c sbb::Container class which is an essential part of the
 * memory management system used by StochBB. Usually the C++ programmer needs to keep track of all
 * objects still in use and is responsibe to free unneeded objects to avoid memory leaks. This can
 * be a difficult taks when dealing with complex structured objects cross referencing eachother.
 * To ease the usage of StochBB, a mark and sweep garbage colector is implemened which keeps track
 * of all objects being directly or indirectly reachable and freeing all unreachable objects. For
 * this memory management system to work, it is necessary to treat all container like values
 * although they represent references to objects allocated on the heap.
 *
 * The central class of StochBB is @c sbb::Var, representing a random variable. This could be a
 * simple random variable having a specified distribution (see @c sbb::AtomicVar) or a random
 * variable that is derived from others like @c sbb::Chain, @c sbb::Minimum, @c sbb::Maximum or
 * @c sbb::Mixture. All random variables have probability density function attached. You can access
 * it using the @c sbb::Var::density method which returns a @c sbb::Density object.
 *
 * All @c sbb::Density objects have two methods, @c sbb::Density::eval evaluating the probability
 * density function and @c sbb::Density::evalCDF evaluating the cumulative density or probability
 * function. Assembling a system of random variables and evaluate their PDFs or CDFs is straight
 * forward. Sampling, however, is not that trivial and is described below in some detail.
 *
 * \section apirv Assembling random processes
 * In a first step, one may define a new gamma-distributed random variable with shape \f$k=10\f$
 * and scale \f$\theta=100\f$ as
 * \code
 *  #include <stochbb/api.h>
 *  using namespace sbb;
 *
 *  // [...]
 *
 *  Var X1 = AtomicVar::gamma(10, 100);
 * \endcode
 *
 * Its PDF can then be evaluated as on a regular grid in \f$[0,1000)\f$ with 1000 grid points with
 * \code
 *  // [...]
 *
 *  Eigen::VectorXd pdf(1000);
 *  X1.density().eval(0, 1000, pdf);
 * \endcode
 * The result of the evaluation is stored into the vector @c pdf. There are only very few basic
 * or atomic random variable types defined in StochBB:
 *
 *   Constructor | Parameters | Process description
 *   --- | --- | ---
 *   @c sbb::AtomicVar::delta | delay | A constant delay or a process with a fixed waiting time.
 *   @c sbb::AtomicVar::unif | a, b | A process with a uniform-distributed waiting time.
 *   @c sbb::AtomicVar::norm | mu, sigma | A process with a normal-distributed waiting.
 *   @c sbb::AtomicVar::gamma | k, theta | A process with a gamma-distributed waiting time.
 *
 * More complex processes can be derived by combining atomic random variables.
 *
 * \subsection apichain Sums of random variables
 * The most basic derived random variable is a @c sbb::Chain. This type represents the "chaining" of
 * random processes. For example, given two random variables \f$X_1\f$ and \f$X_2\f$ and the random
 * process \f$X_1\f$ should trigger the next process \f$X_2\f$, the resulting random process \f$Y\f$
 * can be desribed by the simple sum of the processes \f$X_1\f$ and \f$X_2\f$ as \f$Y=X_1+X_2\f$.
 * That is, the time needed to complete both processes sequentially is simply the sum of the times
 * needed for each process. Such a chain can be constructed using the overloaded @c + operator
 * or the @c sbb::chain function. For example
 * \code
 *  #include <stochbb/api.hh>
 *  using namespace sbb;
 *
 *  // [...]
 *
 *  Var X1 = AtomicVar::gamma(10,100);
 *  Var X2 = AtomicVar::gamma(20, 50);
 *  Var Y = X1 + X2;
 * \endcode
 *
 * \subsection apiminmax Minimum and Maximum of random variables
 * Another simple derived random variable is the @c sbb::Maximum or @c sbb::Minimum class. As the
 * names suggest, they represent the maximum or minimum of a set of random variables. They can be
 * created using the overloaded standard library function @c std::min and @c std::max or the
 * @c sbb::minimum and @c sbb::maximum functions. The latter take a vector of random variables.
 * \code
 *  #include <stochbb/api.hh>
 *  using namespace sbb;
 *
 *  // [...]
 *
 *  Var X1 = AtomicVar::gamma(10,100);
 *  Var X2 = AtomicVar::gamma(20, 50);
 *  Var Y = std::max(X1, X2);
 * \endcode
 *
 * \subsection apimix Mixtures of random variables
 * Similar to the @c sbb::Minimum or @c sbb::Maximum a mixtrue of random variables can be
 * constructed using the @c sbb::mix function. This function takes a vector of weights and a vector
 * of corresponding random variables. Such a mixture can be considered as a random process which
 * randomly selects the outcome of a set of other random processes, where the probability of
 * selecting a specific process is given by the weight assigned to each process.
 * \code
 *  #include <stochbb/api.hh>
 *  using namespace sbb;
 *
 *  // [...]
 *
 *  Var X1 = AtomicVar::gamma(10,100);
 *  Var X2 = AtomicVar::gamma(20, 50);
 *  // Assemble vector of variables
 *  std::vector<Var> vars; var.push_back(X1), vars.push_back(X2);
 *  // Assemble vector of weights
 *  std::vector<double> weights; var.push_back(1), vars.push_back(2);
 *  // Construct mixture
 *  Var Y = sbb::mix(weights, vars);
 * \endcode
 * In the example above, the random process \f$Y\f$ will select the outcome of \f$X_1\f$ with
 * a probability of \f$\frac{1}{3}\f$ and the outcome of \f$X_2\f$ with probability
 * \f$\frac{2}{3}\f$.
 *
 * \subsection apicomp Compound random variables
 * An important class of derived random processes are compound proceses. There the parameters of the
 * distribution of a random variable are themselfs random variables. That is
 * \f[
 *  X \sim f(x|A)\,,\quad A \sim g(a|\theta)\,,
 * \f]
 * where the random variable \f$X\f$ is distributed as \f$f(x|A)\f$, parametrized by \f$A\f$,
 * where \f$A\f$ itself is a random variable distributed as \f$g(a|\theta)\f$, parametrized by
 * \f$\theta\f$. Compound random variables are created using factory methods provided by the
 * @c sbb::Compound class. For example, the @c sbb::Compound::norm factory method constructs
 * a compount-normal distributed random variable, where both the mean and the standard deviation
 * can by any other random variable. For example
 * \code
 *  *  #include <stochbb/api.hh>
 *  using namespace sbb;
 *
 *  // [...]
 *
 *  Var mu = AtomicVar::gamma(10,100);
 *  Var sigma = AtomicVar::delta(10);
 *  Var cnorm = Compound::norm(mu, sigma);
 * \endcode
 * instantiates a compound-normal distributed random variable, where the mean is gamma-distributed
 * while the standard deviation is fixed (implemented by a delta distribution).
 *
 *
 * \section apisample Sample random variables
 * As mentioned above, sampling from a complex process efficiently is not trivial. First of all, is
 * must be ensured that all atomic random variables are samples only once. Otherwise, two dependent
 * random variables may be sampled as independent. Moreover, the results of derived random variables
 * should be chached for efficiency. StochBB provides a separate class that implements a proper
 * sampler for a system of random processes, the @c sbb::ExactSampler class. This class allows to
 * sample from several possibly depdendent random variables simultaneously. Upon construction, the
 * set of random variables to sample is specified. A sample from these random variables can then be
 * obtained by the @c sbb::ExactSampler::sample method.
 * \code
 *  #include <stochbb/api.hh>
 *  using namespace sbb;
 *
 *  // [...]
 *
 *  Var X1 = AtomicVar::gamma(10,100);
 *  Var X2 = AtomicVar::gamma(20, 50);
 *  Var Y = std::min(X1, X2);
 *  // Assemble vector of variables
 *  std::vector<Var> vars; vars.push_back(X1), vars.push_back(X2); vars.push_back(Y);
 *  // Construct sampler
 *  ExactSampler sampler(vars);
 *  // Get 1000 samples
 *  Eigen::MatrixXd samples(3, 1000);
 *  sampler.sample(samples);
 * \endcode
 * The @c sbb::ExactSampler::sample method takes a single @c Eigen::Matrix where each column
 * represents the random variable given to the constructor and each row an independent sample.
 *
 * For very large systems, sampling may become pretty slow. Particularily if one is only interested
 * in the marignal distribution of single random variables. For these cases a approximative sampler
 * for single random variables is provided, the @c sbb::MarginalSampler. This sampler uses an
 * approximation of the inverse of the cummulative distribution function of a random variable
 * to draw samples.
 * \code
 *  *  #include <stochbb/api.hh>
 *  using namespace sbb;
 *
 *  // [...]
 *
 *  Var X1 = AtomicVar::gamma(10,100);
 *  Var X2 = AtomicVar::gamma(20, 50);
 *  Var Y = std::min(X1, X2);
 *  // Sample from Y on [0,500] in 1000 steps
 *  MarginalSampler sampler(Y, 0, 500, 1000);
 *  // Get 1000 samples
 *  Eigen::VectorXd samples(1000);
 *  sampler.sample(samples);
 * \endcode
 * In this example, a @c MarginalSampler is constructed for the random variable Y. Using an
 * approximation of its CDF on the interval \f$[0,500)\f$ using 1000 steps. Then, the
 * @c sbb::MarginalSampler::sample method is used to obtain a sample.
 */

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
class MixtureObj;
class CompoundObj;
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
  /** Packs the given objects. Takes the reference of the object. */
  explicit Container(Object *obj);
  //! @endcond

public:
  /** Empty constructor. */
  Container();

  //! @cond internal
  /** Copy constructor. */
  Container(const Container &other);
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
  static AtomicVar delta(double delay, const std::string &name="");
  /** Constructs a uniformly distributed random variable. */
  static AtomicVar unif(double a, double b, const std::string &name="");
  /** Constructs a normal distributed random variable. */
  static AtomicVar norm(double mu, double sigma, const std::string &name="");
  /** Constructs a gamma distributed random variable. */
  static AtomicVar gamma(double k, double theta, const std::string &name="");

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


/** This class implements a mixture of random variables.
 * A mixture is a random process that selects one of its children with a certain probability. */
class Mixture: public DerivedVar
{
public:
  //! @cond internal
  /** Object type of the container.*/
  typedef MixtureObj ObjectType;

public:
  /** Packs the given @c MixtureObj instance. */
  Mixture(MixtureObj *obj);
  //! @endcond

  /** Constructs a mixture from the given weights and variables. */
  Mixture(const std::vector<double> &weights, const std::vector<Var> &variables, const std::string &name="");

  //! @cond internal
  /** Copy constructor. */
  Mixture(const Mixture &other);
  /** Assignment operator. */
  Mixture &operator =(const Mixture &other);
  //! @endcond

  /** Retruns the weight of the i-th variable. */
  double weight(size_t i) const;

protected:
  //! @cond internal
  /** The reference to the @c MixtureObj instance. */
  MixtureObj *_mixture;
  //! @endcond
};


/** Represents a compound of random variables.
 * Compounds of random variables are formed by treating the parameters of the distribution of
 * one random variable as random too. For example, assume that the random variable
 *  \f$X|A \sim f(x|A)\f$ and \f$A \sim f(a|\theta)\f$. Then, the distribution of \f$X\f$ is
 * \f$f(x|\theta) = \int f(x|a)\,f(a|\theta)\,da\f$. */
class Compound: public DerivedVar
{
public:
  //! @cond internal
  /** The object type of the container. */
  typedef CompoundObj ObjectType;

public:
  /** Packs the given CompoundObj instance. */
  Compound(CompoundObj *obj);
  /** Copy constructor. */
  Compound(const Compound &other);
  /** Assignment operator. */
  Compound &operator =(const Compound &other);
  //! @endcond

  /** Constructs a compound random variable from a normal distribution where
   * \f$\mu\f$ and \f$\sigma\f$ are given by the specified random variables. */
  static Compound norm(const Var &mu, const Var &sigma, const std::string &name="");
  /** Constructs a compound random variable from a gamma distribution where
   * \f$k\f$ and \f$\theta\f$ are given by the specified random variables. */
  static Compound gamma(const Var &k, const Var &theta, const std::string &name="");

protected:
  //! @cond internal
  /** The reference to the @c CompoundObj instance. */
  CompoundObj *_compound;
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

  /** Evaluates the PDF of the selected output variables and stores the results into the given
   * matrix. The matrix gets resized. Each column represents a output variable where the first
   * column is time. Each row represetns a time-point at which the PDFs are evaluated. */
  void evalPDF(Eigen::MatrixXd &out) const;
  /** Evaluates the CDF of the selected output variables and stores the results into the given
   * matrix. The matrix gets resized. Each column represents a output variable where the first
   * column is time. Each row represetns a time-point at which the CDFs are evaluated. */
  void evalCDF(Eigen::MatrixXd &out) const;
  /** Samples from the selected output variables and stores the results into the given
   * matrix. The matrix gets resized. Each column represents a output variable and each row
   * represetns a sample of the variables. */
  void sample(Eigen::MatrixXd &out) const;

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

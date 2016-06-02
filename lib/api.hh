/** @defgroup api Application programming interface
 * Here you find all classes and functions that are part of the public API.
 *
 * StochBB implements a simple mark & sweep garbage collector (GC, @ref mem) for freeing all
 * unreachable objects. For this GC to work, all directly reachable objects are held within
 * containers (classes derived from @c stochbb::Container). These container classes form the API of
 * StochBB and are listed below. */

#ifndef __SBB_API_HH__
#define __SBB_API_HH__

#include "exception.hh"
#include <vector>
#include <Eigen/Eigen>
#include "object.hh"


namespace stochbb {

// Forward declarations
class DensityObj;
class DistributionObj;
class DeltaDistributionObj;
class UniformDistributionObj;
class NormalDistributionObj;
class GammaDistributionObj;
class InvGammaDistributionObj;
class WeibullDistributionObj;
class StudtDistributionObj;
class AtomicDensityObj;
class VarObj;
class AtomicVarObj;
class DerivedVarObj;
class AffineTrafoObj;
class ChainObj;
class MaximumObj;
class MinimumObj;
class MixtureObj;
class ConditionalObj;
class CondChainObj;
class CompoundObj;
class ExactSamplerObj;
class MarginalSamplerObj;


/** Base class of all container classes.
 * @ingroup api */
class Container
{
public:
  /** The object type of the container. */
  typedef Object ObjectType;

protected:
  /** Packs the given objects. Takes the reference of the object. */
  explicit Container(Object *obj);

public:
  /** Empty constructor. */
  Container();

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

  /** Returns @c true if the container is empty. */
  bool isNull() const;
  /** Returns @c true if the container is not empty. */
  bool isValid() const;

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
  /** Empty constructor. */
  Density();
  /** Packs the given @c DensityObj and takes the reference. */
  Density(DensityObj *obj);
  /** Copy constructor. */
  Density(const Density &other);
  /** Assigment operator. */
  Density &operator=(const Density &other);

  /** Evaluates the density (PDF) on a regular grid \f$[Tmin, Tmax)\f$ where the number
   * of grid points is specified via the length of the output vector @c out. The results are
   * stored into the output vector. */
  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;

  /** Evaluates the probability function (CDF) on a regular grid \f$[Tmin, Tmax)\f$ where the number
   * of grid points is specified via the length of the output vector @c out. The results are
   * stored into the output vector. */
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;

  /** Returns an affine transform of the density. */
  Density affine(double scale, double shift) const;

  /** Returns an estimate of the \f$\alpha\f$-quantiles. */
  void rangeEst(double alpha, double &a, double &b) const;

  /** Provides a lexicographical order for all densities. */
  int compare(const Density &other) const;

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


/** Base class of all distributions (families of PDFs).
 * @ingroup api */
class Distribution: public Container
{
public:
  /** The object type of the container. */
  typedef DistributionObj ObjectType;

public:
  /** Packs the given @c DistributionObj and takes the reference. */
  Distribution(DistributionObj *obj);
  /** Copy constructor. */
  Distribution(const Distribution &other);
  /** Assigment operator. */
  Distribution &operator=(const Distribution &other);

  /** Retruns a weak reference to the @c DensityObj. */
  inline DistributionObj *operator *() const { return _distribution; }
  /** Retruns a weak reference to the @c DensityObj. */
  inline DistributionObj *operator ->() const { return _distribution; }
  /** Retruns a new reference to the @c DensityObj. */
  inline DistributionObj *ref() const { _object->ref(); return _distribution; }

  /** Returns the number of parameters. */
  size_t nParams() const;
  /** Evaluates the probability density function with the specified parameters on a regular grid
   * in \f$[T_{min}, T_{max})\f$ and stores the result into @c out. */
  void pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  /** Evaluates the probability function with the specified parameters on a regular grid
   * in \f$[T_{min}, T_{max})\f$ and stores the result into @c out. */
  void cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
           const Eigen::Ref<const Eigen::VectorXd> params) const;
  /** Returns the quantiles (@c lower, @c upper) for the given probability @c p with the specified
   * parameters. */
  void quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const;
  /** Changes the given set of parameters such that the distribution is an affine transformed. */
  void affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const;
  /** Draws some samples from the distribution with the specified parameters. */
  void sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const;
  /** Draws some samples from the distribution with the specified parameters. */
  double sample(const Eigen::Ref<const Eigen::VectorXd> params) const;

protected:
  /** The reference to the @c DistributionObj instance. */
  DistributionObj *_distribution;
};


/** Represents the family of delta distributions.
 * @ingroup api */
class DeltaDistribution: public Distribution
{
public:
  /** The object type of the container. */
  typedef DeltaDistributionObj ObjectType;

public:
  /** Packs the given @c DeltaDistributionObj and takes the reference. */
  DeltaDistribution(DeltaDistributionObj *obj);
  /** Copy constructor. */
  DeltaDistribution(const DeltaDistribution &other);
  /** Assigment operator. */
  DeltaDistribution &operator=(const DeltaDistribution &other);

protected:
  /** The reference to the @c DeltaDistributionObj instance. */
  DeltaDistributionObj *_delta;
};


/** Represents the family of uniform distributions.
 * @ingroup api */
class UniformDistribution: public Distribution
{
public:
  /** The object type of the container. */
  typedef UniformDistributionObj ObjectType;

public:
  /** Packs the given @c UniformDistributionObj and takes the reference. */
  UniformDistribution(UniformDistributionObj *obj);
  /** Copy constructor. */
  UniformDistribution(const UniformDistribution &other);
  /** Assigment operator. */
  UniformDistribution &operator=(const UniformDistribution &other);

protected:
  /** The reference to the @c UniformDistributionObj instance. */
  UniformDistributionObj *_uniform;
};


/** Represents the family of normal distributions.
 * @ingroup api */
class NormalDistribution: public Distribution
{
public:
  /** The object type of the container. */
  typedef NormalDistributionObj ObjectType;

public:
  /** Packs the given @c NormalDistributionObj and takes the reference. */
  NormalDistribution(NormalDistributionObj *obj);
  /** Copy constructor. */
  NormalDistribution(const NormalDistribution &other);
  /** Assigment operator. */
  NormalDistribution &operator=(const NormalDistribution &other);

protected:
  /** The reference to the @c NormalDistributionObj instance. */
  NormalDistributionObj *_normal;
};


/** Represents the family of gamma distributions.
 * @ingroup api */
class GammaDistribution: public Distribution
{
public:
  /** The object type of the container. */
  typedef GammaDistributionObj ObjectType;

public:
  /** Packs the given @c GammaDistributionObj and takes the reference. */
  GammaDistribution(GammaDistributionObj *obj);
  /** Copy constructor. */
  GammaDistribution(const GammaDistribution &other);
  /** Assigment operator. */
  GammaDistribution &operator=(const GammaDistribution &other);

protected:
  /** The reference to the @c GammaDistributionObj instance. */
  GammaDistributionObj *_gamma;
};


/** Represents the family of inverse gamma distributions.
 * @ingroup api */
class InvGammaDistribution: public Distribution
{
public:
  /** The object type of the container. */
  typedef InvGammaDistributionObj ObjectType;

public:
  /** Packs the given @c InvGammaDistributionObj and takes the reference. */
  InvGammaDistribution(InvGammaDistributionObj *obj);
  /** Copy constructor. */
  InvGammaDistribution(const InvGammaDistribution &other);
  /** Assigment operator. */
  InvGammaDistribution &operator=(const InvGammaDistribution &other);

protected:
  /** The reference to the @c InvGammaDistributionObj instance. */
  InvGammaDistributionObj *_invgamma;
};


/** Represents the family of Weibull distributions.
 * @ingroup api */
class WeibullDistribution: public Distribution
{
public:
  /** The object type of the container. */
  typedef WeibullDistributionObj ObjectType;

public:
  /** Packs the given @c WeibullDistributionObj and takes the reference. */
  WeibullDistribution(WeibullDistributionObj *obj);
  /** Copy constructor. */
  WeibullDistribution(const WeibullDistribution &other);
  /** Assigment operator. */
  WeibullDistribution &operator=(const WeibullDistribution &other);

protected:
  /** The reference to the @c WeibullDistributionObj instance. */
  WeibullDistributionObj *_weibull;
};


/** Represents the family of Student's t distributions.
 * @ingroup api */
class StudtDistribution: public Distribution
{
public:
  /** The object type of the container. */
  typedef StudtDistributionObj ObjectType;

public:
  /** Packs the given @c StudtDistributionObj and takes the reference. */
  StudtDistribution(StudtDistributionObj *studt);
  /** Copy constructor. */
  StudtDistribution(const StudtDistribution &other);
  /** Assigment operator. */
  StudtDistribution &operator=(const StudtDistribution &other);

protected:
  /** The reference to the @c StudtDistributionObj instance. */
  StudtDistributionObj *_studt;
};


/** The densities of @c AtomicVar.
 * This is the base class of all densities of atomic random variables, that is a random variable
 * that does not depend on others. */
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

  /** Returns the distribution object of the density. */
  Distribution distribution() const;
  /** Returns the i-th parameter. */
  double parameter(size_t i) const;

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


/** A generic random variable derived from a chosen density. */
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

protected:
  /** Holds the @c GenericVarObj. */
  AtomicVarObj *_genericRV;
};


/** Base class of all derived random variables.
 * @ingroup api */
class DerivedVar: public Var
{
public:
  /** The object type of the container. */
  typedef DerivedVarObj ObjectType;

public:
  /** Packs the @c DerivedVarObj instance. */
  DerivedVar(DerivedVarObj *obj);
  /** Copy constructor. */
  DerivedVar(const DerivedVar &other);
  /** Assignment operator. */
  DerivedVar &operator =(const DerivedVar &other);

  /** Returns the number of random variables the variable depends on directly. */
  size_t numVariables() const;
  /** Returns the i-th random variable the chain depends on directly. */
  Var variable(size_t idx) const;

protected:
  /** Holds a reference to the @c DerivedVarObj instance. */
  DerivedVarObj *_derived_var;
};


/** Represents an affine transformation of another random variable.
 * That is, \f$Y = a\,X+b\f$.
 * @ingroup api */
class AffineTrafo: public DerivedVar
{
public:
  /** Object type of the container. */
  typedef AffineTrafoObj ObjectType;

public:
  /** Packs the given AffineTrafoObj. */
  AffineTrafo(AffineTrafoObj *obj);
  /** Copy constructor. */
  AffineTrafo(const AffineTrafo &other);
  /** Assignment operator. */
  AffineTrafo &operator =(const AffineTrafo &other);

  /** Returns the scale of the transformation. */
  double scale() const;
  /** Returns the shift of the transformation. */
  double shift() const;

protected:
  /** Hols a reference to the @c AffineTrafoObj instance. */
  AffineTrafoObj *_affine;
};


/** Implements the sum of the independent random variables \f$X_i, i=1,\dots,N\f$,
 *  \f$Y = \sum_i X_i\f$.
 *
 * This class is called @c Chain as it implements a chain of random processes. Assume the
 * independent random processes \f$X_1\f$ and \f$X_2\f$ and \f$X_1\f$ triggers \f$X_2\f$,
 * than the time of the second process (\f$X_2\f$) to complete is simply the sum of times
 * of both processes.
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


/** This class implements a mixture of random variables.
 * A mixture is a random process that selects one of its children with a certain probability.
 * @ingroup api */
class Mixture: public DerivedVar
{
public:
  /** Object type of the container.*/
  typedef MixtureObj ObjectType;

public:
  /** Packs the given @c MixtureObj instance. */
  Mixture(MixtureObj *obj);
  /** Constructs a mixture from the given weights and variables. */
  Mixture(const std::vector<double> &weights, const std::vector<Var> &variables, const std::string &name="");
  /** Copy constructor. */
  Mixture(const Mixture &other);
  /** Assignment operator. */
  Mixture &operator =(const Mixture &other);

  /** Retruns the weight of the i-th variable. */
  double weight(size_t i) const;

protected:
  /** The reference to the @c MixtureObj instance. */
  MixtureObj *_mixture;
};


/** This random variable represents a simple conditional mixture of two random variables.
 * Assume \f$X_1, X_2\f$ and \f$Y_1\f$ as well as \f$X_1, X_2\f$ and \f$Y_2\f$ are mutually
 * indepdent random variables. Then, the random variable \f$Z\f$ defined as
 * \f[
 *  Z = \begin{cases}
 *        Y_1 & if\, X_1<X_2 \\
 *        Y_2 & else,
 *       \end{cases}
 * \f]
 * is implemented by this class.
 * @ingroup api */
class Conditional: public DerivedVar
{
public:
  /** Object type of the container.*/
  typedef ConditionalObj ObjectType;

public:
  /** Constructs a conditional mixture from the given random variables. */
  Conditional(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2, const std::string &name="");
  /** Packs the given @c ConditionalObj instance. */
  Conditional(ConditionalObj *obj);
  /** Copy constructor. */
  Conditional(const Conditional &other);
  /** Assignment operator. */
  Conditional &operator=(const Conditional &other);

protected:
  /** Holds a reference to the @c ConditionalObj instance. */
  ConditionalObj *_conditional;
};


/** This random variable represents a conditional chain of two random variables.
 * Assume \f$X_1, X_2\f$ and \f$Y_1\f$ as well as \f$X_1, X_2\f$ and \f$Y_2\f$ are mutually
 * indepdent random variables. Then, the random variable \f$Z\f$ defined as
 * \f[
 *  Z = \begin{cases}
 *        X_1 + Y_1 & if\, X_1<X_2 \\
 *        X_2 + Y_2 & else,
 *       \end{cases}
 * \f]
 * is implemented by this class.
 * @ingroup api */
class CondChain: public DerivedVar
{
public:
  /** Object type of the container. */
  typedef CondChainObj ObjectType;

public:
  /** Constructs a conditional chain random variable from the given RVs. */
  CondChain(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2, const std::string &name="");
  /** Packs the @c CondChainObj instance. */
  CondChain(CondChainObj *obj);
  /** Copy constructor. */
  CondChain(const CondChain &other);
  /** Assignment operator. */
  CondChain &operator=(const CondChain &other);

protected:
  /** Holds a weak reference to the @c CondChainObj instance. */
  CondChainObj *_condchain;
};


/** Represents a compound of random variables.
 * Compounds of random variables are formed by treating the parameters of the distribution of
 * one random variable as random too. For example, assume that the random variable
 *  \f$X|A \sim f(x|A)\f$ and \f$A \sim f(a|\theta)\f$. Then, the distribution of \f$X\f$ is
 * \f$f(x|\theta) = \int f(x|a)\,f(a|\theta)\,da\f$.
 * @ingroup api */
class Compound: public DerivedVar
{
public:
  /** The object type of the container. */
  typedef CompoundObj ObjectType;

public:
  /** Packs the given CompoundObj instance. */
  Compound(CompoundObj *obj);
  /** Copy constructor. */
  Compound(const Compound &other);
  /** Assignment operator. */
  Compound &operator =(const Compound &other);

protected:
  /** The reference to the @c CompoundObj instance. */
  CompoundObj *_compound;
};


/** Implements a sampler for several possibly dependent random variables.
 * This class should be used if one is interesed in sampling several random variables
 * simultaneously. If only one variable of a very large system of random variables is sampled,
 * consider also the @c MarginalSampler class.
 * @ingroup api */
class ExactSampler: public Container
{
public:
  /** Object type of the container. */
  typedef ExactSamplerObj ObjectType;

public:
  /** Constructs a sampler for the given random variable. */
  ExactSampler(const Var &X);
  /** Constructs a sampler for the given random variables. */
  ExactSampler(const Var &X1, const Var &X2);
  /** Constructs a sampler for the given random variables. */
  ExactSampler(const Var &X1, const Var &X2, const Var &X3);
  /** Constructs a sampler for the given random variables. */
  ExactSampler(const std::vector<Var> &variables);

  /** Copy constructor. */
  ExactSampler(const ExactSampler &other);
  /** Assignment operator. */
  ExactSampler &operator =(const ExactSampler &other);

  /** Samples from the random variables passed to the constructor.
   * The number of row of @c out specify the number of samples to draw. The number of columns
   * must match the number of random variables passed to the constructor, i.e. each colum
   * corresponds to a random variable passed to the constructor. */
  void sample(Eigen::Ref<Eigen::MatrixXd> out) const;

protected:
  /** A reference to the sampler object. */
  ExactSamplerObj *_sampler;
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
 * more than one mutually dependent variables are sampled simultaneously.
 * @ingroup api */
class MarginalSampler: public Container
{
public:
  /** Object type of the container. */
  typedef MarginalSamplerObj ObjectType;

public:
  /** Constructs a sampler for the specified variable.
   * @param var Specifies the random variable to sample from.
   * @param Tmin Specifies the lower-bound of the grid on which the CDF is evaluated.
   * @param Tmax Specifies the upper-bound of the grid on which the CDF is evaluated.
   * @param steps Specifies the number of grid points. */
  MarginalSampler(const Var &var, double Tmin, double Tmax, size_t steps);

  /** Copy constructor. */
  MarginalSampler(const MarginalSampler &other);
  /** Assignment operator. */
  MarginalSampler &operator =(const MarginalSampler &other);

  /** Samples from the marginal. The number of samples is specified by the number of elements
   * in the vector @c out. */
  void sample(Eigen::Ref<Eigen::VectorXd> out) const;

protected:
  /** Holds a reference to the sampler object. */
  MarginalSamplerObj *_sampler;
};

}

#include "operators.hh"

#endif // __SBB_API_HH__

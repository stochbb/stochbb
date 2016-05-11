/** @defgroup red Analytic density reductions
 * @ingroup internal */

#ifndef REDUCTION_H
#define REDUCTION_H

#include "api.hh"

namespace stochbb {

/** Represents a single reduction rule for the convolution of two densities.
 * @ingroup red */
class ConvolutionReductionRule
{
protected:
  /** Hidden constructor. Automatically adds the rule instance to the singleton instance
   * of @c stochbb::ConvolutionReductions. */
  ConvolutionReductionRule();

public:
  /** Destructor. */
  virtual ~ConvolutionReductionRule();

  /** Returns @c true if the rule can be applied. */
  virtual bool test(const Density &a, const Density &b) const = 0;
  /** Applies the rule to the densities if @c test returned @c true. */
  virtual Density apply(const Density &a, const Density &b) const = 0;
};


/** Implements the reduction of convolutions with a delta distribution.
 * That is \f$\delta(x-a)\ast f(x) \rightarrow f(x-a)\f$.
 * @ingroup red */
class DeltaConvolutionRule: public ConvolutionReductionRule
{
protected:
  /** Hidden constructor. */
  DeltaConvolutionRule();

public:
  bool test(const Density &a, const Density &b) const;
  Density apply(const Density &a, const Density &b) const;

protected:
  /** Singleton instance. */
  static DeltaConvolutionRule *_rule;
};


/** Implements the convolution of normal distributions.
 * That is \f$f(x,\mu_1, \sigma_1^2)\ast
 *  f(x; \mu_2, \sigma^2_2) \rightarrow f(x;\mu_1+\mu_2, \sigma_1^2+\sigma_2^2)\f$.
 * @ingroup red */
class NormalConvolutionRule: public ConvolutionReductionRule
{
protected:
  /** Hidden constructor. */
  NormalConvolutionRule();

public:
  bool test(const Density &a, const Density &b) const;
  Density apply(const Density &a, const Density &b) const;

protected:
  /** Singleton instance. */
  static NormalConvolutionRule *_rule;
};


/** Implements the convolution of gamma distributions.
 * That is \f$f(x,k_1, \theta)\ast f(x; k_2, \theta) \rightarrow f(x;k_1+k_2, \theta)\f$.
 * @ingroup red */
class GammaConvolutionRule: public ConvolutionReductionRule
{
protected:
  /** Hidden constructor. */
  GammaConvolutionRule();

public:
  bool test(const Density &a, const Density &b) const;
  Density apply(const Density &a, const Density &b) const;

protected:
  /** Singleton instance. */
  static GammaConvolutionRule *_rule;
};


/** Holds several reduction rules.
 * @ingroup red */
class ConvolutionReductions
{
protected:
  /** Hidden construtor. */
  ConvolutionReductions();

public:
  /** Destructor. */
  virtual ~ConvolutionReductions();

  /** Adds a rule to the reduction. The ownership is transferred to the instance. */
  void add(ConvolutionReductionRule *rule);

  /** Searches for a rule matching the given densities. Returns @c null if no matching rule can be
   * found. */
  ConvolutionReductionRule *find(const Density &a, const Density &b) const;

  /** Returns the reductions. */
  static ConvolutionReductions &get();

protected:
  /** Holds all rules. */
  std::vector<ConvolutionReductionRule *> _rules;
  /** Singleton instance. */
  static ConvolutionReductions *_reductions;
};


/** The base class of all compound distribution reduction rules.
 * @ingroup red */
class CompoundReductionRule
{
protected:
  /** Hidden constructor. The constructor adds the rule automatically to the @c CompoundReductions
   * instance. */
  CompoundReductionRule();

public:
  /** Destructor. */
  virtual ~CompoundReductionRule();

  /** Returns @c true if the rule can be applied on the given density. */
  virtual bool test(const Density &a) const = 0;
  /** Applies the rule on the given density and returns the reduced density if @c test returned
   * @c true. */
  virtual Density apply(const Density &a) const = 0;
};


/** A rule reducing trivial delta-compound densities. That is
 * \f$X \sim f,\, Y\sim\delta(y-X)\rightarrow X\f$.
 * @ingroup red */
class DeltaCompoundRule: public CompoundReductionRule
{
protected:
  /** Hidden constructor. */
  DeltaCompoundRule();

public:
  bool test(const Density &a) const;
  Density apply(const Density &a) const;

protected:
  /** Sinleton instance. The rule gets instantiated statically. */
  static DeltaCompoundRule *_instance;
};


/** Reduces compound normal distributions with some fixed sigma to a convolution. That is
 * \f$\mu \sim f,\, X \sim \mathcal{N}(\mu, \sigma^2)\rightarrow
 * X \sim f\ast \mathcal{N}(0,\sigma^2)\f$.
 * @ingroup red */
class NormalCompoundRule: public CompoundReductionRule
{
protected:
  /** Hidden constructor, the rule gets instantiated statically. */
  NormalCompoundRule();

public:
  bool test(const Density &a) const;
  Density apply(const Density &a) const;

protected:
  /** Singleton instance. */
  static NormalCompoundRule *_instance;
};


/** Holds several reduction rules for compound distributions.
 * @ingroup red */
class CompoundReductions
{
protected:
  /** Hidden constructor, use the factor method @c get to obtain an instance. */
  CompoundReductions();

public:
  /** Destructor. */
  virtual ~CompoundReductions();

  /** Adds a rule. */
  void add(CompoundReductionRule *rule);
  /** Finds a rule matching the given density. If no matching rule can be found, @c null is
   * returned. */
  CompoundReductionRule *find(const Density &a) const;

  /** Gets the singleton instance for the compound reductions. */
  static CompoundReductions &get();

protected:
  /** Holds the list of rules. */
  std::vector<CompoundReductionRule *> _rules;
  /** The singeton instance. */
  static CompoundReductions *_instance;
};

}

#endif // REDUCTION_H

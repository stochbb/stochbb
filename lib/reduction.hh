#ifndef REDUCTION_H
#define REDUCTION_H

#include "api.hh"

namespace stochbb {

/** Represents a single reduction rule for the convolution of two densities. */
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
  /** Applies the rule to the densities. */
  virtual Density apply(const Density &a, const Density &b) const = 0;
};


/** Implements the reduction of convolutions with a delta distribution.
 * That is \f$\delta(x-a)\ast f(x) \rightarrow f(x-a)\f$. */
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
 *  f(x; \mu_2, \sigma^2_2) \rightarrow f(x;\mu_1+\mu_2, \sigma_1^2+\sigma_2^2)\f$. */
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
 * That is \f$f(x,k_1, \theta)\ast f(x; k_2, \theta) \rightarrow f(x;k_1+k_2, \theta)\f$. */
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


/** Holds several reduction rules. */
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


class CompoundReductionRule
{
protected:
  CompoundReductionRule();

public:
  virtual ~CompoundReductionRule();

  virtual bool test(const Density &a) const = 0;
  virtual Density apply(const Density &a) const = 0;
};


class DeltaCompoundRule: public CompoundReductionRule
{
protected:
  DeltaCompoundRule();

public:
  bool test(const Density &a) const;
  Density apply(const Density &a) const;

protected:
  static DeltaCompoundRule *_instance;
};


class NormalCompoundRule: public CompoundReductionRule
{
protected:
  NormalCompoundRule();

public:
  bool test(const Density &a) const;
  Density apply(const Density &a) const;

protected:
  static NormalCompoundRule *_instance;
};


class CompoundReductions
{
protected:
  CompoundReductions();

public:
  virtual ~CompoundReductions();

  void add(CompoundReductionRule *rule);
  CompoundReductionRule *find(const Density &a) const;

  static CompoundReductions &get();

protected:
  std::vector<CompoundReductionRule *> _rules;
  static CompoundReductions *_instance;
};

}

#endif // REDUCTION_H

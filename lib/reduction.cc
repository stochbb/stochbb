#include "reduction.hh"
#include "density.hh"
#include "distribution.hh"
#include "compound.hh"
#include "chain.hh"
#include "logger.hh"
#include "operators.hh"

using namespace stochbb;


/* ******************************************************************************************** *
 * Implementation of DensityReduction
 * ******************************************************************************************** */
ConvolutionReductions *ConvolutionReductions::_reductions = 0;

ConvolutionReductions::ConvolutionReductions()
  : _rules()
{
  // pass...
}

ConvolutionReductions::~ConvolutionReductions() {
  for (size_t i=0; i<_rules.size(); i++) {
    delete _rules[i];
  }
}

void
ConvolutionReductions::add(ConvolutionReductionRule *rule) {
  _rules.push_back(rule);
}

ConvolutionReductionRule *
ConvolutionReductions::find(const Density &a, const Density &b) const {
  for (size_t i=0; i<_rules.size(); i++) {
    if (_rules[i]->test(a,b))
      return _rules[i];
  }
  return 0;
}

ConvolutionReductions &
ConvolutionReductions::get() {
  if (_reductions)
    return *_reductions;
  _reductions = new ConvolutionReductions();
  return *_reductions;
}


/* ******************************************************************************************** *
 * Implementation of ConvolutionReductionRule
 * ******************************************************************************************** */
ConvolutionReductionRule::ConvolutionReductionRule()
{
  // register rule with rulelist
  ConvolutionReductions::get().add(this);
}

ConvolutionReductionRule::~ConvolutionReductionRule() {
  // pass...
}


/* ******************************************************************************************** *
 * Implementation of DeltaConvolutionRule
 * ******************************************************************************************** */
DeltaConvolutionRule *DeltaConvolutionRule::_rule = new DeltaConvolutionRule();

DeltaConvolutionRule::DeltaConvolutionRule()
  : ConvolutionReductionRule()
{
  // pass...
}

bool
DeltaConvolutionRule::test(const Density &a, const Density &b) const {
  AtomicDensityObj *a_atomic = dynamic_cast<AtomicDensityObj *>(*a);
  AtomicDensityObj *b_atomic = dynamic_cast<AtomicDensityObj *>(*b);

  if (a_atomic && dynamic_cast<DeltaDistributionObj *>(*a_atomic->distribution())) {
    return true;
  }

  if (b_atomic && dynamic_cast<DeltaDistributionObj *>(*b_atomic->distribution())) {
    return true;
  }

  return false;
}

Density
DeltaConvolutionRule::apply(const Density &a, const Density &b) const {
  AtomicDensityObj *a_atomic = dynamic_cast<AtomicDensityObj *>(*a);
  AtomicDensityObj *b_atomic = dynamic_cast<AtomicDensityObj *>(*b);

  if (a_atomic && dynamic_cast<DeltaDistributionObj *>(*a_atomic->distribution())) {
    Density res = b.affine(1, a_atomic->parameter(0));
    logDebug() << "Reduce convolution of " << a << " and " << b
               << " to " << res << ".";
    return res;
  }

  Density res = a.affine(1, b_atomic->parameter(0));
  logDebug() << "Reduce convolution of " << a << " and " << b
             << " to " << res << ".";
  return res;
}


/* ******************************************************************************************** *
 * Implementation of NormalConvolutionRule
 * ******************************************************************************************** */
NormalConvolutionRule *NormalConvolutionRule::_rule = new NormalConvolutionRule();

NormalConvolutionRule::NormalConvolutionRule()
  : ConvolutionReductionRule()
{
  // pass...
}

bool
NormalConvolutionRule::test(const Density &a, const Density &b) const {
  AtomicDensityObj *a_atomic = dynamic_cast<AtomicDensityObj *>(*a);
  AtomicDensityObj *b_atomic = dynamic_cast<AtomicDensityObj *>(*b);

  // apply if both densities are normal densities...
  return (a_atomic && dynamic_cast<NormalDistributionObj *>(*a_atomic->distribution()))
      && (b_atomic && dynamic_cast<NormalDistributionObj *>(*b_atomic->distribution()));
}

Density
NormalConvolutionRule::apply(const Density &a, const Density &b) const {
  AtomicDensityObj *a_atomic = dynamic_cast<AtomicDensityObj *>(*a);
  AtomicDensityObj *b_atomic = dynamic_cast<AtomicDensityObj *>(*b);
  double a_mu = a_atomic->parameter(0), b_mu = b_atomic->parameter(0),
      a_sig = a_atomic->parameter(1), b_sig = b_atomic->parameter(1);
  Eigen::VectorXd params(2); params << a_mu+b_mu, std::sqrt(a_sig*a_sig + b_sig*b_sig);
  Density res(new AtomicDensityObj(new NormalDistributionObj(), params));
  logDebug() << "Reduce convolution of " << a << " and " << b
             << " to " << res << ".";
  return res;
}


/* ******************************************************************************************** *
 * Implementation of GammaConvolutionRule
 * ******************************************************************************************** */
GammaConvolutionRule *GammaConvolutionRule::_rule = new GammaConvolutionRule();

GammaConvolutionRule::GammaConvolutionRule()
  : ConvolutionReductionRule()
{
  // pass...
}

bool
GammaConvolutionRule::test(const Density &a, const Density &b) const {
  AtomicDensityObj *a_atomic = dynamic_cast<AtomicDensityObj *>(*a);
  AtomicDensityObj *b_atomic = dynamic_cast<AtomicDensityObj *>(*b);

  // apply if both densities are gamma densities with the same theta...
  return (a_atomic && dynamic_cast<GammaDistributionObj *>(*a_atomic->distribution()))
      && (b_atomic && dynamic_cast<GammaDistributionObj *>(*b_atomic->distribution()))
      && (a_atomic->parameter(1) == b_atomic->parameter(1));
}

Density
GammaConvolutionRule::apply(const Density &a, const Density &b) const {
  AtomicDensityObj *a_atomic = dynamic_cast<AtomicDensityObj *>(*a);
  AtomicDensityObj *b_atomic = dynamic_cast<AtomicDensityObj *>(*b);
  double a_k = a_atomic->parameter(0), b_k = b_atomic->parameter(0),
      a_theta = a_atomic->parameter(1),
      a_shift = a_atomic->parameter(2), b_shift = b_atomic->parameter(2);
  Eigen::VectorXd params(3); params << a_k+b_k, a_theta, a_shift+b_shift;
  Density res(new AtomicDensityObj(new GammaDistributionObj(), params));
  logDebug() << "Reduce convolution of " << a << " and " << b
             << " to " << res << ".";
  return res;
}


/* ******************************************************************************************** *
 * Implementation of CompoundReductions
 * ******************************************************************************************** */
CompoundReductions *CompoundReductions::_instance = 0;

CompoundReductions::CompoundReductions()
  : _rules()
{
  // pass...
}

CompoundReductions::~CompoundReductions() {
  for (size_t i=0; i<_rules.size(); i++)
    delete _rules[i];
  _rules.clear();
}

void
CompoundReductions::add(CompoundReductionRule *rule) {
  _rules.push_back(rule);
}

CompoundReductionRule *
CompoundReductions::find(const Density &a) const {
  for (size_t i=0; i<_rules.size(); i++) {
    if (_rules[i]->test(a))
      return _rules[i];
  }
  return 0;
}

CompoundReductions &
CompoundReductions::get() {
  if (_instance)
    return *_instance;

  _instance = new CompoundReductions();
  return *_instance;
}


/* ******************************************************************************************** *
 * Implementation of CompoundReductionRule
 * ******************************************************************************************** */
CompoundReductionRule::CompoundReductionRule()
{
  CompoundReductions::get().add(this);
}

CompoundReductionRule::~CompoundReductionRule() {
  // pass...
}


/* ******************************************************************************************** *
 * Implementation of DeltaCompoundRule
 * ******************************************************************************************** */
DeltaCompoundRule *DeltaCompoundRule::_instance = new DeltaCompoundRule();

DeltaCompoundRule::DeltaCompoundRule()
  : CompoundReductionRule()
{
  // pass...
}

bool
DeltaCompoundRule::test(const Density &a) const {
  CompoundDensityObj *comp = dynamic_cast<CompoundDensityObj *>(*a);
  return comp && dynamic_cast<DeltaDistributionObj *>(*comp->distribution());
}

Density
DeltaCompoundRule::apply(const Density &a) const {
  CompoundDensityObj *comp = dynamic_cast<CompoundDensityObj *>(*a);
  logDebug() << "Reduce compound " << a
             << " to " << comp->parameter(0) << ".";
  return comp->parameter(0);
}


/* ******************************************************************************************** *
 * Implementation of NormalCompoundRule
 * ******************************************************************************************** */
NormalCompoundRule *NormalCompoundRule::_instance = new NormalCompoundRule();

NormalCompoundRule::NormalCompoundRule()
  : CompoundReductionRule()
{
  // pass...
}

bool
NormalCompoundRule::test(const Density &a) const {
  CompoundDensityObj *comp = dynamic_cast<CompoundDensityObj *>(*a);
  if (! (comp && dynamic_cast<NormalDistributionObj *>(*comp->distribution())))
    return false;
  // if sigma is fixed
  AtomicDensityObj *sigma_atomic = dynamic_cast<AtomicDensityObj *>(*comp->parameter(1));
  return (sigma_atomic && dynamic_cast<DeltaDistributionObj *>(*sigma_atomic->distribution()));
}

Density
NormalCompoundRule::apply(const Density &a) const {
  // given that a is a normal distr. with fixed sigma
  CompoundDensityObj *comp = dynamic_cast<CompoundDensityObj *>(*a);
  // Get mu and sigma densities
  Density mu = comp->parameter(0), sigma = comp->parameter(1);
  AtomicDensityObj *sigma_atomic = dynamic_cast<AtomicDensityObj *>(*sigma);

  // Assemble atomic density N(0,sigma)
  Eigen::VectorXd params(2); params << 0, sigma_atomic->parameter(0);
  Density S(new AtomicDensityObj(new NormalDistributionObj(), params));

  // check if conv. of S & mu can be reduced further
  if (ConvolutionReductionRule *rule = ConvolutionReductions::get().find(S, mu)) {
    Density res = rule->apply(S, mu);
    logDebug() << "Reduce compound " << a << " to " << res << ".";
    return res;
  }

  // If convolution cannot be reduced -> constrcut conv. density
  std::vector<Density> densities; densities.reserve(2);
  densities.push_back(mu);
  densities.push_back(S);

  // done
  Density res(new ConvolutionDensityObj(densities));
  logDebug() << "Reduce compound " << a << " to " << res << ".";
  return res;
}


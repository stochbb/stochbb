#include "api.hh"
#include "density.hh"
#include "distribution.hh"
#include "randomvariable.hh"
#include "affinetrafo.hh"
#include "chain.hh"
#include "minmax.hh"
#include "mixture.hh"
#include "conditional.hh"
#include "compound.hh"
#include "exactsampler.hh"
#include "marginalsampler.hh"


using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of Container
 * ********************************************************************************************* */
Container::Container()
  : _object(0)
{
  // pass...
}

Container::Container(Object *obj)
  : _object(obj)
{
  // pass...
}

Container::Container(const Container &other)
  : _object(other._object)
{
  other._object->ref();
}

Container::~Container() {
  // Decrement reference count
  if (_object)
    _object->unref();
  _object = 0;
}

const Container &
Container::operator =(const Container &other)
{
  Object *old_obj = _object;
  _object = other._object;
  // reference new object
  if (_object) { _object->ref(); }
  // unreference old object:
  if (old_obj) { old_obj->unref(); }
  // done.
  return *this;
}

bool
Container::isNull() const {
  return 0 == _object;
}

bool
Container::isValid() const {
  return 0 != _object;
}


/* ********************************************************************************************* *
 * Implementation of RandomVariable container
 * ********************************************************************************************* */
Var::Var()
  : Container(), _randomVariable(0)
{
  // pass...
}

Var::Var(VarObj *obj)
  : Container(obj), _randomVariable(obj)
{
  // pass...
}

Var::Var(const Var &other)
  : Container(other), _randomVariable(other._randomVariable)
{
  // pass...
}

Var &
Var::operator =(const Var &other) {
  Container::operator =(other);
  _randomVariable = other._randomVariable;
  return *this;
}

Density
Var::density() const {
  return _randomVariable->density();
}

bool
Var::dependsOn(const Var &other) const {
  return _randomVariable->dependsOn(other);
}

bool
Var::mutuallyIndep(const Var &other) const {
  return _randomVariable->mutuallyIndep(other);
}

const std::string &
Var::name() const {
  return _randomVariable->name();
}

void
Var::setName(const std::string &name) {
  _randomVariable->setName(name);
}


/* ********************************************************************************************* *
 * Implementation of GenericRandomVariable container
 * ********************************************************************************************* */
AtomicVar::AtomicVar(AtomicVarObj *obj)
  : Var(obj), _genericRV(obj)
{
  // pass...
}

AtomicVar::AtomicVar(const AtomicDensity &density, const std::string &name)
  : Var(new AtomicVarObj(*density, name)),
    _genericRV(static_cast<AtomicVarObj *>(_randomVariable))
{
  // pass...
}

AtomicVar::AtomicVar(const AtomicVar &other)
  : Var(other), _genericRV(other._genericRV)
{
  // pass...
}

AtomicVar &
AtomicVar::operator =(const AtomicVar &other) {
  Var::operator =(other);
  _genericRV = other._genericRV;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of Density container
 * ********************************************************************************************* */
Density::Density()
  : Container(), _density(0)
{
  // pass...
}

Density::Density(DensityObj *obj)
  : Container(obj), _density(obj)
{
  // pass...
}

Density::Density(const Density &other)
  : Container(other), _density(other._density)
{
  // pass...
}

Density &
Density::operator =(const Density &other) {
  Container::operator =(other);
  _density = other._density;
  return *this;
}

void
Density::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  _density->eval(Tmin, Tmax, out);
}

void
Density::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  _density->evalCDF(Tmin, Tmax, out);
}

Density
Density::affine(double scale, double shift) const {
  return _density->affine(scale, shift);
}

void
Density::rangeEst(double alpha, double &a, double &b) const {
  _density->rangeEst(alpha, a, b);
}

int
Density::compare(const Density &other) const {
  if (_density)
    return _density->compare(**other);
  return -1;
}


/* ********************************************************************************************* *
 * Implementation of Distribution container
 * ********************************************************************************************* */
Distribution::Distribution(DistributionObj *obj)
  : Container(obj), _distribution(obj)
{
  // pass...
}

Distribution::Distribution(const Distribution &other)
  : Container(other), _distribution(other._distribution)
{
  // pass...
}

Distribution &
Distribution::operator=(const Distribution &other) {
  Container::operator=(other);
  _distribution = other._distribution;
  return *this;
}

size_t
Distribution::nParams() const {
  return _distribution->nParams();
}

void
Distribution::pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  return _distribution->pdf(Tmin, Tmax, out, params);
}

void
Distribution::cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  return _distribution->cdf(Tmin, Tmax, out, params);
}

void
Distribution::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const{
  _distribution->quantile(lower, upper, p, params);
}

void
Distribution::affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const {
  return _distribution->affine(scale, shift, params);
}

void
Distribution::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  _distribution->sample(out, params);
}

double
Distribution::sample(const Eigen::Ref<const Eigen::VectorXd> params) const {
  Eigen::VectorXd out(1);
  _distribution->sample(out, params);
  return out(0);
}


/* ********************************************************************************************* *
 * Implementation of DeltaDistribution container
 * ********************************************************************************************* */
DeltaDistribution::DeltaDistribution(DeltaDistributionObj *delta)
  : Distribution(delta), _delta(delta)
{
  // pass...
}

DeltaDistribution::DeltaDistribution(const DeltaDistribution &other)
  : Distribution(other), _delta(other._delta)
{
  // pass...
}

DeltaDistribution &
DeltaDistribution::operator= (const DeltaDistribution &other)
{
  Distribution::operator=(other);
  _delta = other._delta;
  return *this;
}

/* ********************************************************************************************* *
 * Implementation of UniformDistribution container
 * ********************************************************************************************* */
UniformDistribution::UniformDistribution(UniformDistributionObj *uniform)
  : Distribution(uniform), _uniform(uniform)
{
  // pass...
}

UniformDistribution::UniformDistribution(const UniformDistribution &other)
  : Distribution(other), _uniform(other._uniform)
{
  // pass...
}

UniformDistribution &
UniformDistribution::operator= (const UniformDistribution &other)
{
  Distribution::operator=(other);
  _uniform = other._uniform;
  return *this;
}

/* ********************************************************************************************* *
 * Implementation of NormalDistribution container
 * ********************************************************************************************* */
NormalDistribution::NormalDistribution(NormalDistributionObj *normal)
  : Distribution(normal), _normal(normal)
{
  // pass...
}

NormalDistribution::NormalDistribution(const NormalDistribution &other)
  : Distribution(other), _normal(other._normal)
{
  // pass...
}

NormalDistribution &
NormalDistribution::operator= (const NormalDistribution &other)
{
  Distribution::operator=(other);
  _normal = other._normal;
  return *this;
}

/* ********************************************************************************************* *
 * Implementation of GammaDistribution container
 * ********************************************************************************************* */
GammaDistribution::GammaDistribution(GammaDistributionObj *gamma)
  : Distribution(gamma), _gamma(gamma)
{
  // pass...
}

GammaDistribution::GammaDistribution(const GammaDistribution &other)
  : Distribution(other), _gamma(other._gamma)
{
  // pass...
}

GammaDistribution &
GammaDistribution::operator= (const GammaDistribution &other)
{
  Distribution::operator=(other);
  _gamma = other._gamma;
  return *this;
}

/* ********************************************************************************************* *
 * Implementation of InvGammaDistribution container
 * ********************************************************************************************* */
InvGammaDistribution::InvGammaDistribution(InvGammaDistributionObj *invgamma)
  : Distribution(invgamma), _invgamma(invgamma)
{
  // pass...
}

InvGammaDistribution::InvGammaDistribution(const InvGammaDistribution &other)
  : Distribution(other), _invgamma(other._invgamma)
{
  // pass...
}

InvGammaDistribution &
InvGammaDistribution::operator= (const InvGammaDistribution &other)
{
  Distribution::operator=(other);
  _invgamma = other._invgamma;
  return *this;
}

/* ********************************************************************************************* *
 * Implementation of WeibullDistribution container
 * ********************************************************************************************* */
WeibullDistribution::WeibullDistribution(WeibullDistributionObj *weibull)
  : Distribution(weibull), _weibull(weibull)
{
  // pass...
}

WeibullDistribution::WeibullDistribution(const WeibullDistribution &other)
  : Distribution(other), _weibull(other._weibull)
{
  // pass...
}

WeibullDistribution &
WeibullDistribution::operator= (const WeibullDistribution &other)
{
  Distribution::operator=(other);
  _weibull = other._weibull;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of StudtDistribution container
 * ********************************************************************************************* */
StudtDistribution::StudtDistribution(StudtDistributionObj *studt)
  : Distribution(studt), _studt(studt)
{
  // pass...
}

StudtDistribution::StudtDistribution(const StudtDistribution &other)
  : Distribution(other), _studt(other._studt)
{
  // pass...
}

StudtDistribution &
StudtDistribution::operator= (const StudtDistribution &other)
{
  Distribution::operator=(other);
  _studt = other._studt;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of AtomicDensity container
 * ********************************************************************************************* */
AtomicDensity::AtomicDensity(AtomicDensityObj *obj)
  : Density(obj), _atomic_density(obj)
{
  // pass...
}

AtomicDensity::AtomicDensity(const AtomicDensity &other)
  : Density(other), _atomic_density(other._atomic_density)
{
  // pass...
}

AtomicDensity &
AtomicDensity::operator =(const AtomicDensity &other) {
  Density::operator =(other);
  _atomic_density = other._atomic_density;
  return *this;
}

void
AtomicDensity::sample(Eigen::VectorXd &out) const {
  _atomic_density->sample(out);
}

Distribution
AtomicDensity::distribution() const {
  return _atomic_density->distribution();
}

double
AtomicDensity::parameter(size_t i) const {
  return _atomic_density->parameter(i);
}


/* ********************************************************************************************* *
 * Implementation of DerivedVar container
 * ********************************************************************************************* */
DerivedVar::DerivedVar(DerivedVarObj *obj)
  : Var(obj), _derived_var(obj)
{
  // pass...
}

DerivedVar::DerivedVar(const DerivedVar &other)
  : Var(other), _derived_var(other._derived_var)
{
  // pass...
}

DerivedVar &
DerivedVar::operator =(const DerivedVar &other) {
  Var::operator =(other);
  _derived_var = other._derived_var;
  return *this;
}

size_t
DerivedVar::numVariables() const {
  return _derived_var->numVariables();
}

Var
DerivedVar::variable(size_t idx) const {
  return _derived_var->variable(idx);
}


/* ********************************************************************************************* *
 * Implementation of AffineTrafo container
 * ********************************************************************************************* */
AffineTrafo::AffineTrafo(AffineTrafoObj *obj)
  : DerivedVar(obj), _affine(obj)
{
  // pass...
}

AffineTrafo::AffineTrafo(const AffineTrafo &other)
  : DerivedVar(other), _affine(other._affine)
{
  // pass...
}

AffineTrafo &
AffineTrafo::operator =(const AffineTrafo &other) {
  DerivedVar::operator =(other);
  _affine = other._affine;
  return *this;
}

double
AffineTrafo::scale() const {
  return _affine->scale();
}

double
AffineTrafo::shift() const {
  return _affine->shift();
}


/* ********************************************************************************************* *
 * Implementation of Chain container
 * ********************************************************************************************* */
Chain::Chain(ChainObj *obj)
  : DerivedVar(obj), _chain(obj)
{
  // pass...
}

Chain::Chain(const std::vector<Var> &variables, const std::string &name)
  : DerivedVar(new ChainObj(variables, name)), _chain(static_cast<ChainObj *>(_randomVariable))
{
  // pass...
}

Chain &
Chain::operator=(const Chain &other) {
  Var::operator=(other);
  _chain = other._chain;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of Maximum container
 * ********************************************************************************************* */
Maximum::Maximum(MaximumObj *obj)
  : DerivedVar(obj), _maximum(obj)
{
  // pass...
}

Maximum::Maximum(const std::vector<Var> &variables, const std::string &name)
  : DerivedVar(new MaximumObj(variables, name)), _maximum(static_cast<MaximumObj *>(_randomVariable))
{
  // pass...
}

Maximum::Maximum(const Maximum &other)
  : DerivedVar(other), _maximum(other._maximum)
{
  // pass...
}

Maximum &
Maximum::operator =(const Maximum &other) {
  Var::operator =(other);
  _maximum = other._maximum;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of Minimum container
 * ********************************************************************************************* */
Minimum::Minimum(MinimumObj *obj)
  : DerivedVar(obj), _minimum(obj)
{
  // pass...
}

Minimum::Minimum(const std::vector<Var> &variables)
  : DerivedVar(new MinimumObj(variables)), _minimum(static_cast<MinimumObj *>(_randomVariable))
{
  // pass...
}

Minimum::Minimum(const Minimum &other)
  : DerivedVar(other), _minimum(other._minimum)
{
  // pass...
}

Minimum &
Minimum::operator =(const Minimum &other) {
  Var::operator =(other);
  _minimum = other._minimum;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of Mixture container
 * ********************************************************************************************* */
Mixture::Mixture(MixtureObj *obj)
  : DerivedVar(obj), _mixture(obj)
{
  // pass...
}

Mixture::Mixture(const std::vector<double> &weights, const std::vector<Var> &variables, const std::string &name)
  : DerivedVar(new MixtureObj(weights, variables, name)), _mixture(static_cast<MixtureObj *>(_object))
{
  // pass...
}

Mixture::Mixture(const Mixture &other)
  : DerivedVar(other), _mixture(other._mixture)
{
  // pass...
}

Mixture &
Mixture::operator =(const Mixture &other) {
  DerivedVar::operator =(other);
  _mixture = other._mixture;
  return *this;
}

double
Mixture::weight(size_t i) const {
  return _mixture->weight(i);
}


/* ********************************************************************************************* *
 * Implementation of Conditional container
 * ********************************************************************************************* */
Conditional::Conditional(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2, const std::string &name)
  : DerivedVar(new ConditionalObj(X1, X2, Y1, Y2, name)),
    _conditional(static_cast<ConditionalObj *>(_object))
{
  // pass...
}

Conditional::Conditional(ConditionalObj *obj)
  : DerivedVar(obj), _conditional(obj)
{
  // pass...
}

Conditional::Conditional(const Conditional &other)
  : DerivedVar(other), _conditional(other._conditional)
{
  // pass...
}

Conditional &
Conditional::operator =(const Conditional &other) {
  DerivedVar::operator =(other);
  _conditional = other._conditional;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of CondChain container
 * ********************************************************************************************* */
CondChain::CondChain(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2, const std::string &name)
  : DerivedVar(new CondChainObj(X1, X2, Y1, Y2, name)),
    _condchain(static_cast<CondChainObj *>(_object))
{
  // pass...
}

CondChain::CondChain(CondChainObj *obj)
  : DerivedVar(obj), _condchain(obj)
{
  // pass...
}

CondChain::CondChain(const CondChain &other)
  : DerivedVar(other), _condchain(other._condchain)
{
  // pass...
}

CondChain &
CondChain::operator =(const CondChain &other) {
  DerivedVar::operator =(other);
  _condchain = other._condchain;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of Compound container
 * ********************************************************************************************* */
Compound::Compound(CompoundObj *obj)
  : DerivedVar(obj), _compound(obj)
{
  // pass...
}

Compound::Compound(const Compound &other)
  : DerivedVar(other), _compound(other._compound)
{
  // pass...
}

Compound &
Compound::operator =(const Compound &other) {
  DerivedVar::operator =(other);
  _compound = other._compound;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of ExactSampler container
 * ********************************************************************************************* */
ExactSampler::ExactSampler(const Var &X)
  : Container(new ExactSamplerObj(X)), _sampler(static_cast<ExactSamplerObj *>(_object))
{
  // pass...
}

ExactSampler::ExactSampler(const Var &X1, const Var &X2)
  : Container(new ExactSamplerObj(X1, X2)), _sampler(static_cast<ExactSamplerObj *>(_object))
{
  // pass...
}

ExactSampler::ExactSampler(const Var &X1, const Var &X2, const Var &X3)
  : Container(new ExactSamplerObj(X1, X2, X3)), _sampler(static_cast<ExactSamplerObj *>(_object))
{
  // pass...
}

ExactSampler::ExactSampler(const std::vector<Var> &variables)
  : Container(new ExactSamplerObj(variables)), _sampler(static_cast<ExactSamplerObj *>(_object))
{
  // pass...
}

ExactSampler::ExactSampler(const ExactSampler &other)
  : Container(other), _sampler(other._sampler)
{
  // pass...
}

ExactSampler &
ExactSampler::operator =(const ExactSampler &other) {
  Container::operator =(other);
  _sampler = other._sampler;
  return *this;
}

void
ExactSampler::sample(Eigen::Ref<Eigen::MatrixXd> out) const {
  _sampler->sample(out);
}


/* ********************************************************************************************* *
 * Implementation of MarginalSampler container
 * ********************************************************************************************* */
MarginalSampler::MarginalSampler(const Var &var, double Tmin, double Tmax, size_t steps)
  : Container(new MarginalSamplerObj(var, Tmin, Tmax, steps)),
    _sampler(static_cast<MarginalSamplerObj *>(_object))
{
  // pass...
}

MarginalSampler::MarginalSampler(const MarginalSampler &other)
  : Container(other), _sampler(other._sampler)
{
  // pass...
}

MarginalSampler &
MarginalSampler::operator =(const MarginalSampler &other) {
  Container::operator =(other);
  _sampler = other._sampler;
  return *this;
}

void
MarginalSampler::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  _sampler->sample(out);
}

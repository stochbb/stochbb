#include "api.hh"
#include "xmlparser.hh"
#include "density.hh"
#include "randomvariable.hh"
#include "affinetrafo.hh"
#include "chain.hh"
#include "minmax.hh"
#include "mixture.hh"
#include "compound.hh"
#include "simulation.hh"
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
  if (0 != _object) { _object->unref(); }
  _object = 0;
}

const Container &
Container::operator =(const Container &other)
{
  Object *old_obj = _object;
  _object = other._object;
  // reference new object
  if (0 != _object) { _object->ref(); }
  // unreference old object:
  if (0 != old_obj) { old_obj->unref(); }
  // done.
  return *this;
}

bool
Container::isNull() const {
  return 0 == _object;
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

AtomicVar
AtomicVar::delta(double delay, const std::string &name) {
  return AtomicVarObj::delta(delay, name);
}

AtomicVar
AtomicVar::unif(double a, double b, const std::string &name) {
  return AtomicVarObj::unif(a,b, name);
}

AtomicVar
AtomicVar::norm(double mu, double sigma, const std::string &name) {
  return AtomicVarObj::norm(mu, sigma, name);
}

AtomicVar
AtomicVar::gamma(double k, double theta, const std::string &name) {
  return AtomicVarObj::gamma(k, theta, name);
}


/* ********************************************************************************************* *
 * Implementation of Density container
 * ********************************************************************************************* */
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

Compound
Compound::norm(const Var &mu, const Var &sigma, const std::string &name) {
  return CompoundObj::norm(mu, sigma, name);
}

Compound
Compound::gamma(const Var &k, const Var &theta, const std::string &name) {
  return CompoundObj::gamma(k, theta, name);
}


/* ********************************************************************************************* *
 * Implementation of Simulation container
 * ********************************************************************************************* */
Simulation::Simulation()
  : Container(new SimulationObj()), _simulation(static_cast<SimulationObj *>(_object))
{
  // pass...
}

Simulation::Simulation(SimulationObj *object)
  : Container(object), _simulation(object)
{
  // pass...
}

Simulation::Simulation(const Simulation &other)
  : Container(other), _simulation(other._simulation)
{
  // pass...
}

Simulation &
Simulation::operator =(const Simulation &other) {
  Container::operator=(other);
  _simulation = other._simulation;
  return *this;
}

Simulation
Simulation::fromXml(const std::string &filename) {
  XmlParser parser;
  return parser.parse(filename.c_str());
}

bool
Simulation::hasVar(const std::string &id) const {
  return _simulation->hasVar(id);
}

Var
Simulation::var(const std::string &id) const {
  return _simulation->var(id);
}

void
Simulation::addVar(const std::string &id, Var &var) const {
  _simulation->addVar(id, *var);
}

double
Simulation::tMin() const {
  return _simulation->tMin();
}

void
Simulation::setTMin(double tMin) {
  _simulation->setTMin(tMin);
}

double
Simulation::tMax() const {
  return _simulation->tMax();
}

void
Simulation::setTMax(double tMax) {
  _simulation->setTMax(tMax);
}

size_t
Simulation::steps() const {
  return _simulation->steps();
}

void
Simulation::setSteps(size_t steps) const {
  _simulation->setSteps(steps);
}

size_t
Simulation::numOutputVars() const {
  return _simulation->numOutputVars();
}

Var
Simulation::outputVar(size_t idx) const {
  return _simulation->outputVar(idx);
}


void
Simulation::addOutputVar(const Var &var) {
  return _simulation->addOutputVar(var);
}


void
Simulation::evalPDF(Eigen::MatrixXd &out) const {
  return _simulation->evalPDF(out);
}

void
Simulation::evalCDF(Eigen::MatrixXd &out) const {
  return _simulation->evalCDF(out);
}

void
Simulation::sample(Eigen::MatrixXd &out) const {
  return _simulation->sample(out);
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

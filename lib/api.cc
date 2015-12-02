#include "api.hh"
#include "xmlparser.hh"
#include "density.hh"
#include "randomvariable.hh"
#include "chain.hh"
#include "minmax.hh"


using namespace sbb;


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
  // pass...
}

Container::~Container() {
  // Decrement reference count
  if (0 != _object) { _object->unref(); }
  _object = 0;
}

const Container &
Container::operator =(const Container &other)
{
  // Unregister current object:
  if (0 != _object) { _object->unref(); }
  _object = other._object;
  // Register new object
  if (0 != _object) { _object->ref(); }
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
GenericVar::GenericVar(GenericVarObj *obj)
  : Var(obj), _genericRV(obj)
{
  // pass...
}

GenericVar::GenericVar(const Density &density, const std::string &name)
  : Var(new GenericVarObj(*density, name)),
    _genericRV(static_cast<GenericVarObj *>(_randomVariable))
{
  // pass...
}

GenericVar::GenericVar(const GenericVar &other)
  : Var(other), _genericRV(other._genericRV)
{
  // pass...
}

GenericVar &
GenericVar::operator =(const GenericVar &other) {
  Var::operator =(other);
  _genericRV = other._genericRV;
  return *this;
}

GenericVar
GenericVar::delta(double delay) {
  return GenericVarObj::delta(delay);
}

GenericVar
GenericVar::unif(double a, double b) {
  return GenericVarObj::unif(a,b);
}

GenericVar
GenericVar::norm(double mu, double sigma) {
  return GenericVarObj::norm(mu, sigma);
}

GenericVar
GenericVar::gamma(double k, double theta) {
  return GenericVarObj::gamma(k, theta);
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
Density::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  _density->eval(Tmin, Tmax, out);
}

void
Density::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  _density->evalCDF(Tmin, Tmax, out);
}


/* ********************************************************************************************* *
 * Implementation of Chain container
 * ********************************************************************************************* */
Chain::Chain(ChainObj *obj)
  : Var(obj), _chain(obj)
{
  // pass...
}

Chain::Chain(const Var &a, const Var &b, const std::string &name)
  : Var(new ChainObj(*a, *b, name)), _chain(static_cast<ChainObj *>(_randomVariable))
{
  // pass...
}

Chain &
Chain::operator=(const Chain &other) {
  Var::operator=(other);
  _chain = other._chain;
  return *this;
}

size_t
Chain::numVariables() const {
  return _chain->numVariables();
}

Var
Chain::variable(size_t idx) const {
  return _chain->variable(idx);
}


/* ********************************************************************************************* *
 * Implementation of Maximum container
 * ********************************************************************************************* */
Maximum::Maximum(MaximumObj *obj)
  : Var(obj), _maximum(obj)
{
  // pass...
}

Maximum::Maximum(const Var &a, const Var &b, const std::string &name)
  : Var(new MaximumObj(*a, *b, name)), _maximum(static_cast<MaximumObj *>(_randomVariable))
{
  // pass...
}

Maximum::Maximum(const Maximum &other)
  : Var(other), _maximum(other._maximum)
{
  // pass...
}

Maximum &
Maximum::operator =(const Maximum &other) {
  Var::operator =(other);
  _maximum = other._maximum;
  return *this;
}

size_t
Maximum::numVariables() const {
  return _maximum->numVariables();
}

Var
Maximum::variable(size_t idx) const {
  return _maximum->variable(idx);
}


/* ********************************************************************************************* *
 * Implementation of Minimum container
 * ********************************************************************************************* */
Minimum::Minimum(MinimumObj *obj)
  : Var(obj), _minimum(obj)
{
  // pass...
}

Minimum::Minimum(const Var &a, const Var &b)
  : Var(new MinimumObj(*a, *b)), _minimum(static_cast<MinimumObj *>(_randomVariable))
{
  // pass...
}

Minimum::Minimum(const Minimum &other)
  : Var(other), _minimum(other._minimum)
{
  // pass...
}

Minimum &
Minimum::operator =(const Minimum &other) {
  Var::operator =(other);
  _minimum = other._minimum;
  return *this;
}

size_t
Minimum::numVariables() const {
  return _minimum->numVariables();
}

Var
Minimum::variable(size_t idx) const {
  return _minimum->variable(idx);
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
  _simulation->addVar(id, var);
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

Var outputVar(size_t idx) const {
  Simulation::return _simulation->outputVar(idx);
}


void
Simulation::addOutputVar(const Var &var) {
  return _simulation->addOutputVar(var);
}


void
Simulation::run(Eigen::MatrixXd &out) const {
  return _simulation->run(out);
}


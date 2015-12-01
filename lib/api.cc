#include "api.hh"
#include "xmlparser.hh"


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
  // Register object with GC
  if (0 != _object) {
    GC::get().box(_object);
  }
}

Container::Container(const Container &other)
  : _object(other._object)
{
  // Register object with GC
  if (0 != _object) {
    GC::get().box(_object);
  }
}

Container::~Container() {
  // Unregister object with GC
  if (0 != _object) {
    GC::get().unbox(_object);
  }
  _object = 0;
}

const Container &
Container::operator =(const Container &other)
{
  // Unregister current object:
  if (0 != _object) { GC::get().unbox(_object); }
  _object = other._object;
  // Register new object
  if (0 != _object) { GC::get().box(_object); }
  // done.
  return *this;
}

bool
Container::isNull() const {
  return 0 == _object;
}

void
Container::box(Object *obj) {
  if (0 != _object) { GC::get().unbox(_object); }
  _object = obj;
  if (0 != _object) { GC::get().box(_object); }
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

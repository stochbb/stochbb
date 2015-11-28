#include "randomvariable.hh"

using namespace sbb;



/* ********************************************************************************************* *
 * Implementation of RandomVariableObj
 * ********************************************************************************************* */
RandomVariableObj::RandomVariableObj()
  : Object()
{
  // pass...
}

RandomVariableObj::~RandomVariableObj() {
  // pass...
}

void
RandomVariableObj::mark() {
  if (isMarked()) { return; }
  Object::mark();
  std::set<RandomVariableObj *>::iterator item = _dependencies.begin();
  for (; item != _dependencies.end(); item++) {
    (*item)->mark();
  }
}


/* ********************************************************************************************* *
 * Implementation of GenericRandomVariableObj
 * ********************************************************************************************* */
GenericRandomVariableObj::GenericRandomVariableObj(DensityObj *density)
  : RandomVariableObj(), _density(density)
{
  // pass...
}

GenericRandomVariableObj::~GenericRandomVariableObj() {
  // pass...
}

void
GenericRandomVariableObj::mark() {
  if (isMarked()) { return; }
  RandomVariableObj::mark();
  _density->mark();
}

DensityObj *
GenericRandomVariableObj::density() {
  return _density;
}

GenericRandomVariableObj *
GenericRandomVariableObj::delta(double delay) {
  return new GenericRandomVariableObj(new DeltaDensityObj(delay));
}

GenericRandomVariableObj *
GenericRandomVariableObj::unif(double a, double b) {
  return new GenericRandomVariableObj(new UniformDensityObj(a, b));
}

GenericRandomVariableObj *
GenericRandomVariableObj::norm(double mu, double sigma) {
  return new GenericRandomVariableObj(new NormalDensityObj(mu, sigma));
}


/* ********************************************************************************************* *
 * Implementation of RandomVariable container
 * ********************************************************************************************* */
RandomVariable::RandomVariable(RandomVariableObj *obj)
  : Container(obj), _randomVariable(obj)
{
  // pass...
}

RandomVariable::RandomVariable(const RandomVariable &other)
  : Container(other), _randomVariable(other._randomVariable)
{
  // pass...
}

RandomVariable &
RandomVariable::operator =(const RandomVariable &other) {
  Container::operator =(other);
  _randomVariable = other._randomVariable;
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of GenericRandomVariable container
 * ********************************************************************************************* */
GenericRandomVariable::GenericRandomVariable(GenericRandomVariableObj *obj)
  : RandomVariable(obj), _genericRV(obj)
{
  // pass...
}

GenericRandomVariable::GenericRandomVariable(const Density &density)
  : RandomVariable(new GenericRandomVariableObj(*density)),
    _genericRV(static_cast<GenericRandomVariableObj *>(_randomVariable))
{
  // pass...
}

GenericRandomVariable::GenericRandomVariable(const GenericRandomVariable &other)
  : RandomVariable(other), _genericRV(other._genericRV)
{
  // pass...
}

GenericRandomVariable &
GenericRandomVariable::operator =(const GenericRandomVariable &other) {
  RandomVariable::operator =(other);
  _genericRV = other._genericRV;
  return *this;
}




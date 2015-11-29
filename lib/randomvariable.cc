#include "randomvariable.hh"

using namespace sbb;



/* ********************************************************************************************* *
 * Implementation of RandomVariableObj
 * ********************************************************************************************* */
VarObj::VarObj()
  : Object()
{
  // pass...
}

VarObj::~VarObj() {
  // pass...
}

void
VarObj::mark() {
  if (isMarked()) { return; }
  Object::mark();
  std::set<VarObj *>::iterator item = _dependencies.begin();
  for (; item != _dependencies.end(); item++) {
    (*item)->mark();
  }
}


/* ********************************************************************************************* *
 * Implementation of GenericRandomVariableObj
 * ********************************************************************************************* */
GenericVarObj::GenericVarObj(DensityObj *density)
  : VarObj(), _density(density)
{
  // pass...
}

GenericVarObj::~GenericVarObj() {
  // pass...
}

void
GenericVarObj::mark() {
  if (isMarked()) { return; }
  VarObj::mark();
  _density->mark();
}

DensityObj *
GenericVarObj::density() {
  return _density;
}

GenericVarObj *
GenericVarObj::delta(double delay) {
  return new GenericVarObj(new DeltaDensityObj(delay));
}

GenericVarObj *
GenericVarObj::unif(double a, double b) {
  return new GenericVarObj(new UniformDensityObj(a, b));
}

GenericVarObj *
GenericVarObj::norm(double mu, double sigma) {
  return new GenericVarObj(new NormalDensityObj(mu, sigma));
}

GenericVarObj *
GenericVarObj::gamma(double k, double theta) {
  return new GenericVarObj(new GammaDensityObj(k, theta));
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

GenericVar::GenericVar(const Density &density)
  : Var(new GenericVarObj(*density)),
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




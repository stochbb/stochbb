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


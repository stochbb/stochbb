#include "randomvariable.hh"
#include <sstream>

using namespace sbb;



/* ********************************************************************************************* *
 * Implementation of RandomVariableObj
 * ********************************************************************************************* */
std::set<std::string> VarObj::_var_names;

VarObj::VarObj(const std::string &name)
  : Object(), _name(name)
{
  if (0 == _name.size()) {
    size_t no = 0;
    std::stringstream buffer; buffer << "X" << no;
    while (VarObj::_var_names.count(buffer.str())) {
      buffer.str("X"); buffer << ++no;
    }
    _name = buffer.str();
    VarObj::_var_names.insert(_name);
  }
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
GenericVarObj::GenericVarObj(DensityObj *density, const std::string &name)
  : VarObj(name), _density(density)
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
GenericVarObj::delta(double delay, const std::string &name) {
  return new GenericVarObj(new DeltaDensityObj(delay), name);
}

GenericVarObj *
GenericVarObj::unif(double a, double b, const std::string &name) {
  return new GenericVarObj(new UniformDensityObj(a, b), name);
}

GenericVarObj *
GenericVarObj::norm(double mu, double sigma, const std::string &name) {
  return new GenericVarObj(new NormalDensityObj(mu, sigma), name);
}

GenericVarObj *
GenericVarObj::gamma(double k, double theta, const std::string &name) {
  return new GenericVarObj(new GammaDensityObj(k, theta), name);
}


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
AtomicVarObj::AtomicVarObj(AtomicDensityObj *density, const std::string &name)
  : VarObj(name), _density(density)
{
  // pass...
}

AtomicVarObj::~AtomicVarObj() {
  // pass...
}

void
AtomicVarObj::mark() {
  if (isMarked()) { return; }
  VarObj::mark();
  _density->mark();
}

Density AtomicVarObj::density() {
  _density->ref();
  return _density;
}

AtomicVarObj *
AtomicVarObj::delta(double delay, const std::string &name) {
  return new AtomicVarObj(new DeltaDensityObj(delay), name);
}

AtomicVarObj *
AtomicVarObj::unif(double a, double b, const std::string &name) {
  return new AtomicVarObj(new UniformDensityObj(a, b), name);
}

AtomicVarObj *
AtomicVarObj::norm(double mu, double sigma, const std::string &name) {
  return new AtomicVarObj(new NormalDensityObj(mu, sigma), name);
}

AtomicVarObj *
AtomicVarObj::gamma(double k, double theta, const std::string &name) {
  return new AtomicVarObj(new GammaDensityObj(k, theta), name);
}


/* ********************************************************************************************* *
 * Implementation of DerivedVarObj
 * ********************************************************************************************* */
DerivedVarObj::DerivedVarObj(const std::vector<Var> &variables, const std::string &name)
  : VarObj(name), _variables()
{
  // Get & store variable objects
  _variables.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _variables.push_back(*variables[i]);
  }

  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    if (! this->mutuallyIndep(_variables[i])) {
      AssumptionError err;
      err << "Cannot assemble chain variable, arguments are not mutually independent.";
      throw err;
    }
    // Add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add variable itself
    _dependencies.insert(_variables[i]);
  }
}

void
DerivedVarObj::mark() {
  if (isMarked()) { return; }
  VarObj::mark();
  // all dependent variable objects are marked by VarObj::mark() through _dependencies
}


/* ********************************************************************************************* *
 * Implementation of VarSetObj
 * ********************************************************************************************* */
VarSet::VarSet()
  : std::set<Var>(), _vars()
{
  // pass...
}

VarSet::VarSet(const std::set<Var> &variables)
  : std::set<Var>(), _vars(variables)
{
  // pass...
}

VarSet::VarSet(const std::vector<Var> &variables)
  : std::set<Var>(), _vars()
{
  _vars.insert(variables.begin(), variables.end());
}

VarSet::VarSet(const VarSet &other)
  : std::set<Var>(), _vars(other._vars)
{
  // pass...
}

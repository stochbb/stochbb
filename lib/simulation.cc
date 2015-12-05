#include "simulation.hh"
#include "exception.hh"
#include "api.hh"


using namespace sbb;

/* ********************************************************************************************* *
 * Implementation of ContextObj
 * ********************************************************************************************* */
ContextObj::ContextObj(ContextObj *parent)
  : Object(), _parent(parent), _symbols(), _parameters()
{
  // pass...
}

ContextObj::~ContextObj() {
  // pass...
}

void
ContextObj::mark() {
  if (isMarked()) { return; }
  Object::mark();
  std::map<std::string, VarObj *>::iterator item = _symbols.begin();
  for (; item != _symbols.end(); item++) {
    item->second->mark();
  }
}

bool
ContextObj::hasVar(const std::string &id) const {
  if (_symbols.count(id)) { return true; }
  if (_parent) { return _parent->hasVar(id); }
  return false;
}

Var ContextObj::var(const std::string &id) const {
  if (_symbols.count(id)) {
    _symbols.find(id)->second->ref();
    return _symbols.find(id)->second;
  }
  if (_parent) {
    return _parent->var(id);
  }
  return 0;
}

void
ContextObj::addVar(const std::string &id, const Var &var) {
  if (0 != _symbols.count(id)) {
    ParserError err;
    err << "ParserError: Cannot redefine variable '" << id << "'.";
    throw err;
  }
  _symbols[id] = *var;
}

bool
ContextObj::hasParam(const std::string &id) const {
  return ( 0 != _parameters.count(id) );
}

double
ContextObj::param(const std::string &id) const {
  return _parameters.find(id)->second;
}

void
ContextObj::setParam(const std::string &id, double value) {
  _parameters[id] = value;
}


/* ********************************************************************************************* *
 * Implementation of SimulationObj
 * ********************************************************************************************* */
SimulationObj::SimulationObj()
  : ContextObj(0)
{
  // pass...
}

SimulationObj::~SimulationObj() {
  // pass...
}

void
SimulationObj::mark() {
  if (isMarked()) { return; }
  ContextObj::mark();
  for (size_t i=0; i<_outputVariables.size(); i++) {
    _outputVariables[i]->mark();
  }
}

double
SimulationObj::tMin() const {
  return _tMin;
}

void
SimulationObj::setTMin(double Tmin) {
  _tMin = Tmin;
}

double
SimulationObj::tMax() const {
  return _tMax;
}

void
SimulationObj::setTMax(double tMax) {
  _tMax = tMax;
}

size_t
SimulationObj::steps() const {
  return _steps;
}

void
SimulationObj::setSteps(size_t steps) {
  _steps = steps;
}

void
SimulationObj::addOutputVar(const Var &var) {
  _outputVariables.push_back(*var);
}

void
SimulationObj::evalPDF(Eigen::MatrixXd &out) const {
  size_t N = steps();
  size_t M = _outputVariables.size();
  out.resize(N, M+1);

  double Tmin = tMin(), Tmax = tMax(), dt = (Tmax-Tmin)/N;
  double t=Tmin;

  // Get all PDFs and store them into out at column j+1
  for (size_t j=0; j<M; j++) {
    Eigen::VectorXd pdf(N);
    _outputVariables[j]->density()->eval(Tmin, Tmax, pdf);
    out.col(j+1) = pdf;
  }

  // Store time column
  for (size_t i=0; i<N; i++, t+=dt) {
    out(i,0) = t;
  }
}


void
SimulationObj::evalCDF(Eigen::MatrixXd &out) const {
  size_t N = steps();
  size_t M = _outputVariables.size();
  out.resize(N, M+1);

  double Tmin = tMin(), Tmax = tMax(), dt = (Tmax-Tmin)/N;
  double t=Tmin;

  // Get all CDFs and store them into out at column j+1
  for (size_t j=0; j<M; j++) {
    Eigen::VectorXd cdf(N);
    _outputVariables[j]->density()->evalCDF(Tmin, Tmax, cdf);
    out.col(j+1) = cdf;
  }

  // Store time column
  for (size_t i=0; i<N; i++, t+=dt) {
    out(i,0) = t;
  }
}

void
SimulationObj::sample(Eigen::MatrixXd &out) const {
  size_t N = steps();
  size_t M = _outputVariables.size();
  out.resize(N, M);

  std::vector<Var> vars; vars.reserve(M);
  for (size_t i=0; i<M; i++) {
    _outputVariables[i]->ref();
    vars.push_back(_outputVariables[i]);
  }

  ExactSampler sampler(vars);
  sampler.sample(out);
}

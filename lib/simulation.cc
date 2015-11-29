#include "simulation.hh"

using namespace sbb;

/* ********************************************************************************************* *
 * Implementation of SimulationObj
 * ********************************************************************************************* */
SimulationObj::SimulationObj()
  : Object(), _symbols()
{
  // pass...
}

SimulationObj::~SimulationObj() {
  // pass...
}

void
SimulationObj::mark() {
  if (isMarked()) { return; }
  Object::mark();
  std::map<std::string, VarObj *>::iterator item = _symbols.begin();
  for (; item != _symbols.end(); item++) {
    item->second->mark();
  }
  for (size_t i=0; i<_outputVariables.size(); i++) {
    _outputVariables[i]->mark();
  }
}

bool
SimulationObj::hasVar(const std::string &id) {
  return 0 != _symbols.count(id);
}

VarObj *
SimulationObj::var(const std::string &id) {
  return _symbols[id];
}

void
SimulationObj::addVar(const std::string &id, VarObj *var) {
  _symbols[id] = var;
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

const std::vector<VarObj *> &
SimulationObj::outputVars() const {
  return _outputVariables;
}

void
SimulationObj::addOutputVar(VarObj *var) {
  _outputVariables.push_back(var);
}

void
SimulationObj::run(Eigen::MatrixXd &out) const {
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



#ifndef SIMULATION_HH
#define SIMULATION_HH

#include "randomvariable.hh"
#include <map>
#include <vector>


namespace sbb {

class SimulationObj: public Object
{
public:
  SimulationObj();
  virtual ~SimulationObj();

  virtual void mark();

  bool hasVar(const std::string &id);
  VarObj *var(const std::string &id);
  void addVar(const std::string &id, VarObj *variable);

  double tMin() const;
  void setTMin(double tMin);
  double tMax() const;
  void setTMax(double tMax);
  size_t steps() const;
  void setSteps(size_t steps);
  const std::vector<VarObj *> &outputVars() const;
  void addOutputVar(VarObj *var);

  void run(Eigen::MatrixXd &out) const;

protected:
  std::map<std::string, VarObj *> _symbols;
  double _tMin, _tMax;
  size_t _steps;
  std::vector<VarObj *> _outputVariables;
};


class Simulation: public Container
{
public:
  typedef SimulationObj ObjectType;

public:
  Simulation();
  Simulation(SimulationObj *object);
  Simulation(const Simulation &other);

  Simulation &operator= (const Simulation &other);

  inline bool hasVar(const std::string &id) const {
    return _simulation->hasVar(id);
  }
  inline Var var(const std::string &id) const {
    return _simulation->var(id);
  }
  inline void addVar(const std::string &id, Var &var) const {
    _simulation->addVar(id, *var);
  }

  inline double tMin() const { return _simulation->tMin(); }
  inline void setTMin(double tMin) { _simulation->setTMin(tMin); }
  inline double tMax() const { return _simulation->tMax(); }
  inline void setTMax(double tMax) { _simulation->setTMax(tMax); }
  inline size_t steps() const { return _simulation->steps(); }
  inline void setSteps(size_t steps) const { _simulation->setSteps(steps); }

  inline size_t numOutputVars() const { return _simulation->outputVars().size(); }
  inline Var outputVar(size_t idx) const { return _simulation->outputVars()[idx]; }
  inline void addOutputVar(const Var &var) { return _simulation->addOutputVar(*var); }

  inline void run(Eigen::MatrixXd &out) const { return _simulation->run(out); }

protected:
  SimulationObj *_simulation;
};

}

#endif // SIMULATION_HH

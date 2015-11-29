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

}

#endif // SIMULATION_HH

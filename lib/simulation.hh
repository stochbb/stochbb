#ifndef SIMULATION_HH
#define SIMULATION_HH

#include "randomvariable.hh"
#include <map>
#include <vector>


namespace sbb {


class ContextObj: public Object
{
public:
  ContextObj(ContextObj *parent);
  virtual ~ContextObj();
  virtual void mark();

  bool hasVar(const std::string &id) const;
  VarObj *var(const std::string &id) const;
  void addVar(const std::string &id, VarObj *variable);

  bool hasParam(const std::string &id) const;
  double param(const std::string &id) const;
  void setParam(const std::string &id, double value);

protected:
  ContextObj *_parent;
  std::map<std::string, VarObj *> _symbols;
  std::map<std::string, double> _parameters;
};


class SimulationObj: public ContextObj
{
public:
  SimulationObj();
  virtual ~SimulationObj();

  virtual void mark();

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
  double _tMin, _tMax;
  size_t _steps;
  std::vector<VarObj *> _outputVariables;
};

}

#endif // SIMULATION_HH

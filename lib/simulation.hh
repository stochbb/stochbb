#ifndef SIMULATION_HH
#define SIMULATION_HH

#include "api.hh"
#include "randomvariable.hh"
#include <map>
#include <vector>


namespace sbb {

/** A context mapping unique identifier to random variables and parameters.
 * @ingroup internal */
class ContextObj: public Object
{
public:
  /** Constructor from parent context. */
  ContextObj(ContextObj *parent);
  /** Destructor. */
  virtual ~ContextObj();
  virtual void mark();

  /** Returns @c true if the given identifier is assinged to a random variable. */
  bool hasVar(const std::string &id) const;
  /** Returns the random variable specified by the identifier. */
  Var var(const std::string &id) const;
  /** Adds a random variable to the context. */
  void addVar(const std::string &id, const Var &variable);
  /** Returns @c true if the given parameter is specified. */
  bool hasParam(const std::string &id) const;
  /** Returns the value of the specified parameter. */
  double param(const std::string &id) const;
  /** Sets the parameter. */
  void setParam(const std::string &id, double value);

protected:
  /** The parent context. */
  ContextObj *_parent;
  /** The random variable table. */
  std::map<std::string, VarObj *> _symbols;
  /** The parameter table. */
  std::map<std::string, double> _parameters;
};


/** Represents a simulation.
 * @ingroup internal */
class SimulationObj: public ContextObj
{
public:
  /** Empty constructor. */
  SimulationObj();
  /** Destructor. */
  virtual ~SimulationObj();

  virtual void mark();

  /** Returns the start time. */
  double tMin() const;
  /** Sets the start time. */
  void setTMin(double tMin);

  /** Returns the end time. */
  double tMax() const;
  /** Sets the end time. */
  void setTMax(double tMax);

  /** Returns the number of time steps. */
  size_t steps() const;
  /** Sets the number of time steps. */
  void setSteps(size_t steps);

  /** Returns the number of output variables defined in the simulation. */
  inline size_t numOutputVars() const {
    return _outputVariables.size();
  }
  /** Returns the specified output variable. */
  inline Var outputVar(size_t i) const {
    _outputVariables[i]->ref();
    return _outputVariables[i];
  }
  /** Adds an output variable to the simulation. */
  void addOutputVar(const Var &var);

  /** Evaluates the PDF of the selected output variables. */
  void evalPDF(Eigen::MatrixXd &out) const;
  /** Evaluates the CDF of the selected output variables. */
  void evalCDF(Eigen::MatrixXd &out) const;
  /** Samples from the selected output variables. */
  void sample(Eigen::MatrixXd &out) const;

protected:
  /** Start time. */
  double _tMin;
  /** End time. */
  double _tMax;
  /** Time steps. */
  size_t _steps;
  /** The output variables. */
  std::vector<VarObj *> _outputVariables;
};

}

#endif // SIMULATION_HH

#include "exactsampler.hh"
#include "randomvariable.hh"
#include "chain.hh"
#include "minmax.hh"
#include "exception.hh"
#include <algorithm>

using namespace sbb;


ExactSamplerObj::ExactSamplerObj(const std::vector<Var> &variables)
  : Object(), _outvars(), _queue(), _varmap(), _sampler()
{
  _queue.reserve(2*variables.size());
  _outvars.reserve(variables.size());
  // Assemble variable queue
  for (size_t i=0; i<variables.size(); i++) {
    _add_to_queue(*variables[i]);
    _outvars.push_back(*variables[i]);
  }
}

void
ExactSamplerObj::mark() {
  if (isMarked()) { return; }
  Object::mark();
  for (size_t i=0; i<_queue.size(); i++) {
    _queue[i]->mark();
  }
}

void
ExactSamplerObj::sample(Eigen::MatrixXd &out) {
  // Allocate sample cache
  Eigen::MatrixXd samples(out.rows(), _queue.size());
  // Sample all variables in reverse order
  for (int i=(_queue.size()-1); i>=0; i--) {
    _sampler[i](this, _queue[i], samples);
  }
  // Collect selected variables and store them into out
  for (size_t i=0; i<_outvars.size(); i++) {
    out.col(i) = samples.col(_varmap[_outvars[i]]);
  }
}

void
ExactSamplerObj::_add_to_queue(VarObj *var) {
  // Check if variable is already in the queue
  if (0 != _varmap.count(var)) { return; }
  // If not add to queue
  _varmap[var] = _queue.size();
  _queue.push_back(var);
  // Select sampler for the random variable
  _sampler.push_back(_choose_sampler(var));
  // if the variable is a derived variable -> process direct children recusively
  if (DerivedVarObj *derived = dynamic_cast<DerivedVarObj *>(var)) {
    for (size_t i=0; i<derived->numVariables(); i++) {
      _add_to_queue(*derived->variable(i));
    }
  }
}

ExactSamplerObj::sampler_f
ExactSamplerObj::_choose_sampler(VarObj *var) {
  // dispatch by var type
  if (dynamic_cast<AtomicVarObj *>(var)) {
    return ExactSamplerObj::_sample_atomic;
  } else if (dynamic_cast<ChainObj *>(var)) {
    return ExactSamplerObj::_sample_chain;
  } else if (dynamic_cast<MinimumObj *>(var)) {
    return ExactSamplerObj::_sample_minimum;
  } else if (dynamic_cast<MaximumObj *>(var)) {
    return ExactSamplerObj::_sample_minimum;
  }

  TypeError err;
  err << "Cannot create sampler for RV " << var
      << ": Unknown type.";
  throw err;
}

void
ExactSamplerObj::_sample_atomic(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  Eigen::VectorXd tmp(out.rows());
  static_cast<AtomicVarObj *>(var)->sample(tmp);
  out.col(self->_varmap[var]) = tmp;
}

void
ExactSamplerObj::_sample_chain(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  ChainObj *chain = static_cast<ChainObj *>(var);
  size_t var_idx = self->_varmap[var];
  out.col(var_idx).setZero();
  for (size_t i=0; i<chain->numVariables(); i++) {
    size_t idx = self->_varmap[*chain->variable(i)];
    out.col(var_idx) += out.col(idx);
  }
}

void
ExactSamplerObj::_sample_minimum(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  MinimumObj *min_obj = static_cast<MinimumObj *>(var);
  size_t var_idx = self->_varmap[var];
  out.col(var_idx) = out.col(self->_varmap[*min_obj->variable(0)]);
  for (size_t i=1; i<min_obj->numVariables(); i++) {
    size_t idx = self->_varmap[*min_obj->variable(i)];
    for (int j=0; j<out.rows(); j++) {
      out(j, var_idx) = std::min(out(j, var_idx), out(j, idx));
    }
  }
}

void
ExactSamplerObj::_sample_maximum(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  MaximumObj *max_obj = static_cast<MaximumObj *>(var);
  size_t var_idx = self->_varmap[var];
  out.col(var_idx) = out.col(self->_varmap[*max_obj->variable(0)]);
  for (size_t i=1; i<max_obj->numVariables(); i++) {
    size_t idx = self->_varmap[*max_obj->variable(i)];
    for (int j=0; j<out.rows(); j++) {
      out(j, var_idx) = std::max(out(j, var_idx), out(j, idx));
    }
  }
}

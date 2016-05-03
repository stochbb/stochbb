#include "exactsampler.hh"
#include "randomvariable.hh"
#include "affinetrafo.hh"
#include "chain.hh"
#include "minmax.hh"
#include "mixture.hh"
#include "compound.hh"
#include "conditional.hh"
#include "exception.hh"
#include "rng.hh"
#include "logger.hh"

#include <algorithm>
using namespace stochbb;


inline size_t _find_index(double p, size_t a, size_t b, const std::vector<double> &cdf) {
  while (1 < (b-a)) {
    size_t mid = a+(b-a)/2;
    if (p < cdf[mid]) { b = mid; }
    else { a = mid; }
  }
  return b;
}

ExactSamplerObj::ExactSamplerObj(const Var &X) throw (TypeError)
  : Object(), _outvars(), _queue(), _varmap(), _sampler()
{
  _queue.reserve(2);
  _outvars.reserve(1);
  // Assemble variable queue
  _add_to_queue(*X);
  _outvars.push_back(*X);
}

ExactSamplerObj::ExactSamplerObj(const Var &X1, const Var &X2) throw (TypeError)
  : Object(), _outvars(), _queue(), _varmap(), _sampler()
{
  _queue.reserve(4);
  _outvars.reserve(2);
  // Assemble variable queue
  _add_to_queue(*X1);
  _outvars.push_back(*X1);
  _add_to_queue(*X2);
  _outvars.push_back(*X2);
}

ExactSamplerObj::ExactSamplerObj(const Var &X1, const Var &X2, const Var &X3) throw (TypeError)
  : Object(), _outvars(), _queue(), _varmap(), _sampler()
{
  _queue.reserve(6);
  _outvars.reserve(3);
  // Assemble variable queue
  _add_to_queue(*X1);
  _outvars.push_back(*X1);
  _add_to_queue(*X2);
  _outvars.push_back(*X2);
  _add_to_queue(*X3);
  _outvars.push_back(*X3);
}

ExactSamplerObj::ExactSamplerObj(const std::vector<Var> &variables) throw (TypeError)
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
ExactSamplerObj::sample(Eigen::Ref<Eigen::MatrixXd> out) {
  // Allocate sample cache
  Eigen::MatrixXd samples(out.rows(), _queue.size());
  // Sample all variables in reverse order
  for (size_t i=0; i<_queue.size(); i++) {
    _sampler[i](this, _queue[i], samples);
  }
  // Collect selected variables and store them into out
  for (size_t i=0; i<_outvars.size(); i++) {
    out.col(i) = samples.col(_varmap[_outvars[i]]);
  }
}

void
ExactSamplerObj::_add_to_queue(VarObj *var) throw (TypeError) {
  // Check if variable is already in the queue
  if (0 != _varmap.count(var)) { return; }

  // if the variable is a derived variable -> process direct children recusively
  if (DerivedVarObj *derived = dynamic_cast<DerivedVarObj *>(var)) {
    for (size_t i=0; i<derived->numVariables(); i++) {
      _add_to_queue(*derived->variable(i));
    }
  }

  // If add to queue
  _varmap[var] = _queue.size();
  _queue.push_back(var);
  // Select sampler for the random variable
  _sampler.push_back(_choose_sampler(var));
}

ExactSamplerObj::sampler_f
ExactSamplerObj::_choose_sampler(VarObj *var) throw (TypeError) {
  // dispatch by var type
  if (dynamic_cast<AtomicVarObj *>(var)) {
    return ExactSamplerObj::_sample_atomic;
  } else if (dynamic_cast<AffineTrafoObj *>(var)) {
    return ExactSamplerObj::_sample_affine;
  } else if (dynamic_cast<ChainObj *>(var)) {
    return ExactSamplerObj::_sample_chain;
  } else if (dynamic_cast<MinimumObj *>(var)) {
    return ExactSamplerObj::_sample_minimum;
  } else if (dynamic_cast<MaximumObj *>(var)) {
    return ExactSamplerObj::_sample_minimum;
  } else if (dynamic_cast<MixtureObj *>(var)) {
    return ExactSamplerObj::_sample_mixture;
  } else if (dynamic_cast<ConditionalObj *>(var)) {
    return ExactSamplerObj::_sample_conditional;
  } else if (dynamic_cast<CondChainObj *>(var)) {
    return ExactSamplerObj::_sample_condchain;
  } else if (dynamic_cast<CompoundObj *>(var)) {
    return ExactSamplerObj::_sample_compound;
  }

  TypeError err;
  err << "Cannot create sampler for RV ";
  var->print(err);
  err << ": Unknown type.";
  throw err;
}

void
ExactSamplerObj::_sample_atomic(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  Eigen::VectorXd tmp(out.rows());
  static_cast<AtomicVarObj *>(var)->sample(tmp);
  out.col(self->_varmap[var]) = tmp;
}

void
ExactSamplerObj::_sample_affine(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  AffineTrafoObj *affine = static_cast<AffineTrafoObj *>(var);
  // Get samples from untransformed variable and apply affine transform
  out.col(self->_varmap[var])
      = affine->scale()*out.col(self->_varmap[*affine->variable(0)]).array()+affine->shift();
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

void
ExactSamplerObj::_sample_mixture(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  MixtureObj *mix = static_cast<MixtureObj *>(var);

  // compute cummulative distribution of weights
  std::vector<double> cum; cum.reserve(mix->numVariables());
  for (size_t i=0; i<mix->numVariables(); i++) {
    cum[i] = (i>0) ? (cum[i-1]+mix->weight(i)) : mix->weight(i);
  }
  // normalize cdf
  for (size_t i=0; i<mix->numVariables(); i++) {
    cum[i] /= cum[mix->numVariables()-1];
  }

  // get output column
  size_t var_idx = self->_varmap[var];
  for (int i=0; i<out.rows(); i++) {
    // select a variable randomly
    double p = RNG::unif();
    size_t idx = (p < cum[0]) ? 0 : _find_index(p, 0, mix->numVariables()-1, cum);
    // select sample
    out(i, var_idx) = out(i, self->_varmap[*mix->variable(idx)]);
  }
}

void
ExactSamplerObj::_sample_conditional(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  ConditionalObj *cond = static_cast<ConditionalObj *>(var);
  size_t X1 = self->_varmap[*cond->variable(0)];
  size_t X2 = self->_varmap[*cond->variable(1)];
  size_t Y1 = self->_varmap[*cond->variable(2)];
  size_t Y2 = self->_varmap[*cond->variable(3)];
  size_t Z  = self->_varmap[cond];
  for (int i=0; i<out.size(); i++) {
    out(i, Z) = (out(i,X1)<out(i,X2)) ? out(i,Y1) : out(i,Y2);
  }
}

void
ExactSamplerObj::_sample_condchain(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  CondChainObj *cond = static_cast<CondChainObj *>(var);
  size_t X1 = self->_varmap[*cond->variable(0)];
  size_t X2 = self->_varmap[*cond->variable(1)];
  size_t Y1 = self->_varmap[*cond->variable(2)];
  size_t Y2 = self->_varmap[*cond->variable(3)];
  size_t Z  = self->_varmap[cond];
  for (int i=0; i<out.size(); i++) {
    out(i, Z) = (out(i,X1)<out(i,X2)) ? out(i,X1)+out(i,Y1) : out(i,X2)+out(i,Y2);
  }
}

void
ExactSamplerObj::_sample_compound(ExactSamplerObj *self, VarObj *var, Eigen::MatrixXd &out) {
  CompoundObj *comp = static_cast<CompoundObj *>(var);

  std::vector<size_t> p_idxs; p_idxs.reserve(comp->distribution().nParams());
  // Get indices of parameters
  for (size_t i=0; i<p_idxs.size(); i++) {
    p_idxs[i] = self->_varmap[*comp->parameter(i)];
  }
  size_t Z  = self->_varmap[comp];

  Eigen::VectorXd params(p_idxs.size());
  for (int i=0; i<out.size(); i++) {
    // Assemble sample of parameters
    for (size_t j=0; j<p_idxs.size(); j++) {
      params(j) = out(i, p_idxs[j]);
    }
    // get a single sample with sampled parameters:
    out(i, Z) = comp->distribution().sample(params);
  }
}

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
  : Object(), _outvars(), _queue(), _varmap()
{
  _queue.reserve(2);
  _outvars.reserve(1);
  // Assemble variable queue
  _add_to_queue(*X);
  _outvars.push_back(*X);
}

ExactSamplerObj::ExactSamplerObj(const Var &X1, const Var &X2) throw (TypeError)
  : Object(), _outvars(), _queue(), _varmap()
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
  : Object(), _outvars(), _queue(), _varmap()
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
  : Object(), _outvars(), _queue(), _varmap()
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
  for (size_t i=0; i<_outvars.size(); i++) {
    _outvars[i]->mark();
  }
}

void
ExactSamplerObj::sample(Eigen::Ref<Eigen::MatrixXd> out) {
  // Allocate sample cache
  Eigen::MatrixXd samples(out.rows(), _queue.size());
  // allocate index vector
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> indices(_queue.size());

  // Sample all variables in reverse order
  for (size_t i=0; i<_queue.size(); i++) {
    if (AtomicVarObj *avar = dynamic_cast<AtomicVarObj *>(_queue[i])) {
      avar->sample(samples.col(_varmap[_queue[i]]));
    } else if (DerivedVarObj *dvar = dynamic_cast<DerivedVarObj *>(_queue[i])) {
      // Assemble index vector
      for (size_t j=0; j<dvar->numVariables(); j++) {
        indices(j) = _varmap[*dvar->variable(j)];
      }
      // sample
      dvar->sample(_varmap[_queue[i]], indices.head(dvar->numVariables()), samples);
    } else {
      TypeError err;
      err << "Cannot sample from "; _queue[i]->print(err);
      err << ", unknown variable type.";
      throw err;
    }
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
}

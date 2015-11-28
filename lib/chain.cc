#include "chain.hh"
#include <unsupported/Eigen/FFT>

using namespace sbb;


/* ********************************************************************************************* *
 * Implementation of ConvolutionDensityObj
 * ********************************************************************************************* */
ConvolutionDensityObj::ConvolutionDensityObj(const std::vector<RandomVariableObj *> &variables)
  : DensityObj(), _densities()
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(variables[i]->density());
  }
}

ConvolutionDensityObj::~ConvolutionDensityObj() {
  // pass...
}

void
ConvolutionDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->mark();
  }
}

void
ConvolutionDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  Eigen::FFT<double> fft;
  Eigen::VectorXd tmp1(2*out.size());  tmp1.setZero();
  Eigen::VectorXcd tmp2(2*out.size());
  Eigen::VectorXcd prod(2*out.size()); prod.setOnes();
  // Perform FFT convolution
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->eval(Tmin, Tmax, out);
    tmp1.head(out.size()) = out;
    fft.fwd(tmp2, tmp1); prod.array() *= tmp2.array();
  }
  fft.inv(tmp1, prod);
  out = tmp1.head(out.size());
}

void
ConvolutionDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  this->eval(Tmin, Tmax, out); double dt = (Tmax-Tmin)/out.size();
  /// @todo Implement mid-point method.
  for (int i=1; i<out.size(); i++) {
    out[i] = out[i-1] + out[i]*dt;
  }
}

void
ConvolutionDensityObj::sample(Eigen::VectorXd &out) const {
  Eigen::VectorXd tmp(out.size()); out.setZero();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->sample(tmp); out += tmp;
  }
}


/* ********************************************************************************************* *
 * Implementation of ChainObj
 * ********************************************************************************************* */
ChainObj::ChainObj(RandomVariableObj *a, RandomVariableObj *b)
  : RandomVariableObj(), _variables(), _density(0)
{
  if (ChainObj *a_chain = dynamic_cast<ChainObj *>(a)) {
    _variables = a_chain->variables();
  } else {
    _variables.push_back(a);
  }
  if (ChainObj *b_chain = dynamic_cast<ChainObj *>(b)) {
    _variables.reserve(_variables.size()+b_chain->variables().size());
    for (size_t i=0; i<b_chain->variables().size(); i++) {
      _variables.push_back(b_chain->variables()[i]);
    }
  } else {
    _variables.push_back(b);
  }

  // Create density
  _density = new ConvolutionDensityObj(_variables);

  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    // Add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add variable itself
    _dependencies.insert(_variables[i]);
  }
}

ChainObj::ChainObj(const std::vector<RandomVariableObj *> &variables)
  : RandomVariableObj(), _variables(variables)
{
  _density = new ConvolutionDensityObj(_variables);
  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    // Add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add variable itself
    _dependencies.insert(_variables[i]);
  }
}

ChainObj::~ChainObj() {
  // pass...
}

void
ChainObj::mark() {
  if (isMarked()) { return; }
  RandomVariableObj::mark();
  // mark all variables held
  for (size_t i=0; i<_variables.size(); i++) {
    _variables[i]->mark();
  }
  // mark density
  _density->mark();
}

DensityObj *
ChainObj::density() {
  return _density;
}


/* ********************************************************************************************* *
 * Implementation of Chain container
 * ********************************************************************************************* */
Chain::Chain(const RandomVariable &a, const RandomVariable &b)
  : RandomVariable(new ChainObj(*a, *b)), _chain(static_cast<ChainObj *>(_randomVariable))
{
  // pass...
}

Chain &
Chain::operator=(const Chain &other) {
  RandomVariable::operator=(other);
  _chain = other._chain;
  return *this;
}



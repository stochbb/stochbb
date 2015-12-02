#include "chain.hh"
#include "exception.hh"
#include <unsupported/Eigen/FFT>

using namespace sbb;


/* ********************************************************************************************* *
 * Implementation of ConvolutionDensityObj
 * ********************************************************************************************* */
ConvolutionDensityObj::ConvolutionDensityObj(const std::vector<VarObj *> &variables)
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
  double dt = (Tmax-Tmin)/out.size();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->eval(Tmin, Tmax, out);
    tmp1.head(out.size()) = out*dt;
    fft.fwd(tmp2, tmp1); prod.array() *= tmp2.array();
  }
  fft.inv(tmp1, prod);
  out = tmp1.head(out.size())/dt;
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
ChainObj::ChainObj(VarObj *a, VarObj *b, const std::string &name)
  : VarObj(name), _variables(), _density(0)
{
  // Check if these two RVs are mutually independent
  if (! a->mutuallyIndep(b)) {
    AssumptionError err;
    err << "Cannot assemble chain variable, arguments are not mutually independent.";
    throw err;
  }

  // Flatten chain a
  if (ChainObj *a_chain = dynamic_cast<ChainObj *>(a)) {
    // if a is a chain too -> add its variables directly
    _variables.reserve(a_chain->numVariables());
    for (size_t i=0; i<a_chain->numVariables(); i++) {
      _variables.push_back(*a_chain->variable(i));
    }
  } else {
    // otherwise simply add a
    _variables.push_back(a);
  }

  // Flatten chain b
  if (ChainObj *b_chain = dynamic_cast<ChainObj *>(b)) {
    // if b is a chain too -> add its variables directly
    _variables.reserve(_variables.size()+b_chain->numVariables());
    for (size_t i=0; i<b_chain->numVariables(); i++) {
      _variables.push_back(*b_chain->variable(i));
    }
  } else {
    // otherwise
    _variables.push_back(b);
  }

  // Create density
  _density = new ConvolutionDensityObj(_variables);

  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    // add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add explicit dependency
    _dependencies.insert(_variables[i]);
  }
}

ChainObj::ChainObj(const std::vector<VarObj *> &variables, const std::string &name)
  : VarObj(name), _variables(variables)
{
  _density = new ConvolutionDensityObj(_variables);
  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    if (! this->mutuallyIndep(_variables[i])) {
      AssumptionError err;
      err << "Cannot assemble chain variable, arguments are not mutually independent.";
      throw err;
    }
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
  VarObj::mark();
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


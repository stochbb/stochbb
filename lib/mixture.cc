#include "mixture.hh"
#include "exception.hh"

using namespace stochbb;



/* ********************************************************************************************* *
 * Implementation of MixtrueDensityObj
 * ********************************************************************************************* */
MixtureDensityObj::MixtureDensityObj(const std::vector<double> &weights, const std::vector<DensityObj *> &densities,
                                     double scale, double shift)
  : DensityObj(), _weights(weights), _densities(densities), _scale(scale), _shift(shift)
{
  // pass...
}

MixtureDensityObj::MixtureDensityObj(const std::vector<double> &weights, const std::vector<VarObj *> &variables,
                                     double scale, double shift)
  : DensityObj(), _weights(weights), _densities(), _scale(scale), _shift(shift)
{
  // Normalize weights
  double sum = 0;
  for (size_t i=0; i<_weights.size(); i++) { sum += _weights[i]; }
  for (size_t i=0; i<_weights.size(); i++) { _weights[i] /= sum; }
  // Store densities
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(*variables[i]->density());
  }
}

void
MixtureDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->mark();
  }
}

void
MixtureDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // apply affine transform
  Tmin = (Tmin - _shift)/_scale;
  Tmax = (Tmax - _shift)/_scale;

  out.setZero();
  Eigen::VectorXd tmp(out.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->eval(Tmin, Tmax, tmp);
    out += _weights[i]*tmp/_scale;
  }
}

void
MixtureDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // apply affine transform
  Tmin = (Tmin - _shift)/_scale;
  Tmax = (Tmax - _shift)/_scale;

  out.setZero();
  Eigen::VectorXd tmp(out.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->evalCDF(Tmin, Tmax, tmp);
    out += _weights[i]*tmp;
  }
}

Density
MixtureDensityObj::affine(double scale, double shift) const {
  return new MixtureDensityObj(_weights, _densities, scale*_scale, scale*_shift + shift);
}

void
MixtureDensityObj::rangeEst(double alpha, double &a, double &b) const {
  // Returns an estimate of the quantile as the union of all quantiles
  _densities[0]->rangeEst(alpha, a, b);
  for (size_t i=0; i<_densities.size(); i++) {
    double c,d; _densities[i]->rangeEst(alpha, c, d);
    a = std::min(a, c);
    b = std::max(b, d);
  }
}


/* ********************************************************************************************* *
 * Implementation of MixtrueObj
 * ********************************************************************************************* */
MixtureObj::MixtureObj(const std::vector<double> &weights, const std::vector<Var> &variables,
                       const std::string &name)
  : DerivedVarObj(variables, name), _weights(weights), _density(0)
{
  if (_weights.size() != variables.size()) {
    AssumptionError err;
    err << "Cannot create MixtureObj, number of weights must match number of variables.";
    throw err;
  }

  _density = new MixtureDensityObj(_weights, _variables);
  _density->unref();
}

void
MixtureObj::mark() {
  if (isMarked()) { return; }
  DerivedVarObj::mark();
  if (_density) { _density->mark(); }
}

Density
MixtureObj::density() {
  _density->ref();
  return _density;
}

#include "density.hh"
#include "rng.hh"
#include "math.hh"
#include "logger.hh"
#include "distribution.hh"

#include <typeinfo>
#include <typeindex>
#include <cmath>

using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of DensityObj
 * ********************************************************************************************* */
DensityObj::DensityObj()
  : Object()
{
  // pass...
}

DensityObj::~DensityObj() {
  // pass...
}

void
DensityObj::mark() {
  if (isMarked()) { return; }
  Object::mark();
}

int
DensityObj::compare(const DensityObj &other) const {
  // Same object -> equal
  if (this == &other) { return 0; }

  // compare by type
  if (typeid(*this).before(typeid(other))) { return -1; }
  else if (typeid(other).before(typeid(*this))) { return -1; }
  return 0;
}

void
DensityObj::print(std::ostream &stream) const {
  stream << "<Density #" << (void *)this << ">";
}


/* ********************************************************************************************* *
 * Implementation of GenericAtomicDensityObj
 * ********************************************************************************************* */
AtomicDensityObj::AtomicDensityObj(const Distribution &dist, Eigen::Ref<Eigen::VectorXd> params)
  : DensityObj(), _distribution(*dist), _params(params)
{
  assume(size_t(_params.size()) == _distribution->nParams());
  logDebug() << "Construct " << dist
             << " density with parameters [" << params.transpose() << "].";
}

AtomicDensityObj::~AtomicDensityObj() {
  // pass...
}

void
AtomicDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  if (_distribution) { _distribution->mark(); }
}

Distribution
AtomicDensityObj::distribution() const {
  _distribution->ref();
  return _distribution;
}

size_t
AtomicDensityObj::nParams() const {
  return _distribution->nParams();
}

double
AtomicDensityObj::parameter(size_t i) const {
  return _params[i];
}

void
AtomicDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  _distribution->pdf(Tmin, Tmax, out, _params);
}

void
AtomicDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  _distribution->cdf(Tmin, Tmax, out, _params);
}

Density
AtomicDensityObj::affine(double scale, double shift) const {
  // copy parameters
  Eigen::VectorXd params = _params;
  _distribution->affine(scale, shift, params);
  return new AtomicDensityObj(_distribution, params);
}

void
AtomicDensityObj::rangeEst(double alpha, double &a, double &b) const {
  _distribution->quantile(a, b, alpha, _params);
}

void
AtomicDensityObj::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  _distribution->sample(out, _params);
}

int
AtomicDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // Compare atomic densities
  const AtomicDensityObj *ogen = dynamic_cast<const AtomicDensityObj *>(&other);
  // Compare by distribution type
  if (int res = _distribution->compare(*ogen->_distribution))
    return res;
  // compare by parameters
  for (int i=0; i<_params.size(); i++) {
    if (_params[i] < ogen->_params[i])
      return -1;
    else if (_params[i] > ogen->_params[i])
      return 1;
  }
  return 0;
}

void
AtomicDensityObj::print(std::ostream &stream) const {
  stream << "<AtomicDensity distr="; _distribution->print(stream);
  stream << " params=[" << _params.transpose() << "] #" << this << ">";
}



#include "density.hh"
#include "rng.hh"
#include "math.hh"
#include "logger.hh"

using namespace sbb;


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


/* ********************************************************************************************* *
 * Implementation of DeltaDensityObj
 * ********************************************************************************************* */
DeltaDensityObj::DeltaDensityObj(double delay)
  : DensityObj(), _delay(delay)
{
  logDebug() << "Create DeltaDensity with delay=" << _delay << ".";
}

DeltaDensityObj::~DeltaDensityObj() {
  // pass...
}

void
DeltaDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
}

void
DeltaDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  out.setZero();
  if ((_delay<Tmin) || (_delay>Tmax)) { return; }
  double dt = (Tmax-Tmin)/out.size();
  size_t idx = (_delay-Tmin)/dt;
  out[idx] = 1./dt;
}

void
DeltaDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  out.setZero();
  if ((_delay<Tmin) || (_delay>Tmax)) { return; }
  double dt = (Tmax-Tmin)/out.size();
  size_t idx = (_delay-Tmin)/dt;
  out.tail(out.size()-idx).setConstant(1);
}

void
DeltaDensityObj::sample(Eigen::VectorXd &out) const {
  out.setConstant(_delay);
}


/* ********************************************************************************************* *
 * Implementation of UniformDensityObj
 * ********************************************************************************************* */
UniformDensityObj::UniformDensityObj(double a, double b)
  : DensityObj(), _a(a), _b(b)
{
  logDebug() << "Create UniformDensity with a=" << _a << ", b=" << _b << ".";
}

UniformDensityObj::~UniformDensityObj() {
  // pass...
}

void
UniformDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
}

void
UniformDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = ((t >= _a) && (t <= _b)) ? 1./(Tmax-Tmin) : 0.0;
  }
}

void
UniformDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    if (t<_a) { out[i] = 0.0; }
    else if (t>_b) { out[i] = 1.0; }
    else { out[i] = (t-_a)/(_b-_a); }
  }
}

void
UniformDensityObj::sample(Eigen::VectorXd &out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = RNG::unif(_a, _b);
  }
}


/* ********************************************************************************************* *
 * Implementation of NormalDensityObj
 * ********************************************************************************************* */
NormalDensityObj::NormalDensityObj(double mean, double stddev)
  : DensityObj(), _mu(mean), _sigma(stddev)
{
  logDebug() << "Create NormalDensity with mu=" << _mu << ", sigma=" << _sigma << ".";
}

NormalDensityObj::~NormalDensityObj() {
  // pass...
}

void
NormalDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
}

void
NormalDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = std::exp( -(t-_mu)*(t-_mu)/(2*_sigma*_sigma) ) / (_sigma*std::sqrt(2*M_PI));
  }
}

void
NormalDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = 0.5*(1+std::erf((t-_mu)/(_sigma*std::sqrt(2))));
  }
}

void
NormalDensityObj::sample(Eigen::VectorXd &out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = RNG::norm(_mu, _sigma);
  }
}


/* ********************************************************************************************* *
 * Implementation of GammaDensityObj
 * ********************************************************************************************* */
GammaDensityObj::GammaDensityObj(double k, double theta)
  : DensityObj(), _k(k), _theta(theta)
{
  logDebug() << "Create GammaDensity with k=" << _k << ", theta=" << _theta << ".";
}

GammaDensityObj::~GammaDensityObj() {
  // pass...
}

void
GammaDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
}

void
GammaDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    if (t<=0) { out[i] = 0; }
    else { out[i] = std::exp((_k-1)*std::log(t) - t/_theta -std::lgamma(_k) -_k*std::log(_theta)); }
  }
}

void
GammaDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  double c = std::tgamma(_k);
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = sbb::gamma_li(_k, t/_theta) / c;
  }
}

void
GammaDensityObj::sample(Eigen::VectorXd &out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = RNG::gamma(_k, _theta);
  }
}


#include "density.hh"
#include "rng.hh"

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
  // pass...
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
  // pass...
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
    out[i] = ((t >= _a) && (t <= _b)) ? 1.0 : 0.0;
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
  : DensityObj(), _mean(mean), _stddev(stddev)
{
  // pass...
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
    out[i] = std::exp( -(t-_mean)*(t-_mean)/(2*_stddev*_stddev) ) / (_stddev*std::sqrt(2*M_PI));
  }
}

void
NormalDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = 0.5*(1+std::erf((t-_mean)/(_stddev*std::sqrt(2))));
  }
}

void
NormalDensityObj::sample(Eigen::VectorXd &out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = RNG::norm(_mean, _stddev);
  }
}


/* ********************************************************************************************* *
 * Implementation of Density container
 * ********************************************************************************************* */
Density::Density(DensityObj *obj)
  : Container(obj), _density(obj)
{
  // pass...
}

Density::Density(const Density &other)
  : Container(other), _density(other._density)
{
  // pass...
}

Density &
Density::operator =(const Density &other) {
  Container::operator =(other);
  _density = other._density;
  return *this;
}


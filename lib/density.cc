#include "density.hh"
#include "rng.hh"
#include "math.hh"
#include "logger.hh"
#include <typeinfo>
#include <typeindex>


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
  stream << "<DensityObj #" << (void *)this << ">";
}


/* ********************************************************************************************* *
 * Implementation of AtomicDensityObj
 * ********************************************************************************************* */
AtomicDensityObj::AtomicDensityObj()
  : DensityObj()
{
  // pass...
}

AtomicDensityObj::~AtomicDensityObj() {
  // pass...
}

void
AtomicDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
}


/* ********************************************************************************************* *
 * Implementation of DeltaDensityObj
 * ********************************************************************************************* */
DeltaDensityObj::DeltaDensityObj(double delay)
  : AtomicDensityObj(), _delay(delay)
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
DeltaDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Evaluate delta PDF
  out.setZero();
  if ((_delay<Tmin) || (_delay>Tmax)) { return; }
  double dt = (Tmax-Tmin)/out.size();
  size_t idx = (_delay-Tmin)/dt;
  out[idx] = 1./dt;
}

void
DeltaDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Eval CDF
  out.setZero();
  if ((_delay<Tmin) || (_delay>Tmax)) { return; }
  double dt = (Tmax-Tmin)/out.size();
  size_t idx = (_delay-Tmin)/dt;
  out.tail(out.size()-idx).setConstant(1);
}

void
DeltaDensityObj::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  out.setConstant(_delay);
}

Density
DeltaDensityObj::affine(double scale, double shift) const {
  return new DeltaDensityObj(scale*_delay+shift);
}

int
DeltaDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // Compare delta densities
  const DeltaDensityObj *odelta = dynamic_cast<const DeltaDensityObj *>(&other);
  if (_delay < odelta->_delay) { return -1; }
  else if (_delay > odelta->_delay) { return 1; }
  return 0;
}

void
DeltaDensityObj::print(std::ostream &stream) const {
  stream << "<DeltaDensityObj delay=" << _delay << " #" << (void *)this << ">";
}


/* ********************************************************************************************* *
 * Implementation of UniformDensityObj
 * ********************************************************************************************* */
UniformDensityObj::UniformDensityObj(double a, double b)
  : AtomicDensityObj(), _a(a), _b(b)
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
UniformDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = ((t >= _a) && (t <= _b)) ? 1./(Tmax-Tmin) : 0.0;
  }
}

void
UniformDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    if (t<_a) { out[i] = 0.0; }
    else if (t>_b) { out[i] = 1.0; }
    else { out[i] = (t-_a)/(_b-_a); }
  }
}

void
UniformDensityObj::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = RNG::unif(_a, _b);
  }
}

Density
UniformDensityObj::affine(double scale, double shift) const {
  return new UniformDensityObj(scale*_a+shift, scale*_b+shift);
}

int
UniformDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // Compare uniform densities
  const UniformDensityObj *ounif = dynamic_cast<const UniformDensityObj *>(&other);
  if (_a < ounif->_a) { return -1; }
  else if (_a > ounif->_a) { return 1; }
  if (_b < ounif->_b) { return -1; }
  else if (_b > ounif->_b) { return 1; }
  return 0;
}

void
UniformDensityObj::print(std::ostream &stream) const {
  stream << "<UniformDensityObj a=" << _a << ", b=" << _b << " #" << (void *)this << ">";
}


/* ********************************************************************************************* *
 * Implementation of NormalDensityObj
 * ********************************************************************************************* */
NormalDensityObj::NormalDensityObj(double mean, double stddev)
  : AtomicDensityObj(), _mu(mean), _sigma(stddev)
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
NormalDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = std::exp( -(t-_mu)*(t-_mu)/(2*_sigma*_sigma) ) / (_sigma*std::sqrt(2*M_PI));
  }
}

void
NormalDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = 0.5*(1+std::erf((t-_mu)/(_sigma*std::sqrt(2))));
  }
}

void
NormalDensityObj::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = RNG::norm(_mu, _sigma);
  }
}

Density
NormalDensityObj::affine(double scale, double shift) const {
  return new NormalDensityObj(scale*_mu+shift, scale*_sigma);
}

int
NormalDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // Compare uniform densities
  const NormalDensityObj *onorm= dynamic_cast<const NormalDensityObj *>(&other);
  if (_mu < onorm->_mu) { return -1; }
  else if (_mu > onorm->_mu) { return 1; }
  if (_sigma < onorm->_sigma) { return -1; }
  else if (_sigma > onorm->_sigma) { return 1; }
  return 0;
}

void
NormalDensityObj::print(std::ostream &stream) const {
  stream << "<NormalDensityObj mu=" << _mu << ", sigma=" << _sigma
         << " #" << (void *)this << ">";
}


/* ********************************************************************************************* *
 * Implementation of GammaDensityObj
 * ********************************************************************************************* */
GammaDensityObj::GammaDensityObj(double k, double theta, double shift)
  : AtomicDensityObj(), _k(k), _theta(theta), _shift(shift)
{
  logDebug() << "Create GammaDensity with k=" << _k << ", theta=" << _theta
             << ", shift=" << _shift << "." ;
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
GammaDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform on arguments
  Tmin -= _shift; Tmax -= _shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    if (t<=0) { out[i] = 0; }
    else { out[i] = std::exp((_k-1)*std::log(t) - t/_theta -std::lgamma(_k) -_k*std::log(_theta)); }
  }
}

void
GammaDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform on arguments
  Tmin -= _shift; Tmax -= _shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = sbb::gamma_li(_k, t/_theta);
  }
}

void
GammaDensityObj::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = RNG::gamma(_k, _theta)+_shift;
  }
}

Density
GammaDensityObj::affine(double scale, double shift) const {
  return new GammaDensityObj(_k, scale*_theta, scale*_shift+shift);
}

int
GammaDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // Compare uniform densities
  const GammaDensityObj *ogamma = dynamic_cast<const GammaDensityObj *>(&other);
  if (_k < ogamma->_k) { return -1; }
  else if (_k > ogamma->_k) { return 1; }
  if (_theta < ogamma->_theta) { return -1; }
  else if (_theta > ogamma->_theta) { return 1; }
  if (_shift < ogamma->_shift) { return -1; }
  else if (_shift > ogamma->_shift) { return 1; }
  return 0;
}

void
GammaDensityObj::print(std::ostream &stream) const {
  stream << "<GammaDensityObj k=" << _k << ", theta=" << _theta
         << ", shift=" << _shift << " #" << (void *)this << ">";
}


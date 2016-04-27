#include "density.hh"
#include "rng.hh"
#include "math.hh"
#include "logger.hh"
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
  logDebug() << "Create DeltaDensity with delay=" << _delay << " #" << this << ".";
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

void
DeltaDensityObj::rangeEst(double alpha, double &a, double &b) const {
  a = b = _delay;
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
UniformDensityObj::UniformDensityObj(double a, double b) throw (Error)
  : AtomicDensityObj(), _a(a), _b(b)
{
  logDebug() << "Create UniformDensity with a=" << _a << ", b=" << _b << " #" << this << ".";
  assume(_a<_b);
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
    out[i] = ((t >= _a) && (t <= _b)) ? 1./(_b-_a) : 0.0;
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

void
UniformDensityObj::rangeEst(double alpha, double &a, double &b) const {
  double d = (_b-_a)*alpha/2;
  a = _a+d; b = _b-d;
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
NormalDensityObj::NormalDensityObj(double mean, double stddev) throw (Error)
  : AtomicDensityObj(), _mu(mean), _sigma(stddev)
{
  logDebug() << "Create NormalDensity with mu=" << _mu << ", sigma=" << _sigma
             << " #" << this << ".";
  assume(_sigma > 0);
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
    out[i] = stochbb::dnorm((t-_mu)/_sigma)/_sigma;
  }
}

void
NormalDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::pnorm((t-_mu)/_sigma);
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

void
NormalDensityObj::rangeEst(double alpha, double &a, double &b) const {
  double d = stochbb::qnorm(alpha/2);
  if (alpha<0.5) {
    a = _mu + d*_sigma;
    b = _mu - d*_sigma;
  } else {
    a = _mu - d*_sigma;
    b = _mu + d*_sigma;
  }
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
GammaDensityObj::GammaDensityObj(double k, double theta, double shift) throw (Error)
  : AtomicDensityObj(), _k(k), _theta(theta), _shift(shift)
{
  logDebug() << "Create GammaDensity with k=" << _k << ", theta=" << _theta
             << ", shift=" << _shift << " #" << this << "." ;
  assume(_k > 0);
  assume(_theta > 0);
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
    out[i] = stochbb::dgamma(t, _k, _theta);
  }
}

void
GammaDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform on arguments
  Tmin -= _shift; Tmax -= _shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::pgamma(t, _k, _theta);
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

void
GammaDensityObj::rangeEst(double alpha, double &a, double &b) const {
  a = stochbb::qgamma(alpha/2, _k, _theta);
  b = stochbb::qgamma(1-alpha/2, _k, _theta);
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
  stream << "<GammaDensityObj k=" << _k << ", theta=" << _theta;
  if (_shift) { stream << ", shift=" << _shift; }
  stream << " #" << (void *)this << ">";
}


/* ********************************************************************************************* *
 * Implementation of InvGammaDensityObj
 * ********************************************************************************************* */
InvGammaDensityObj::InvGammaDensityObj(double alpha, double beta, double shift) throw (Error)
  : AtomicDensityObj(), _alpha(alpha), _beta(beta), _shift(shift)
{
  logDebug() << "Create InvGammaDensity with alpha=" << _alpha << ", beta=" << _beta
             << ", shift=" << _shift << " #" << this << "." ;
  assume(_alpha > 0);
  assume(_beta  > 0);
}

InvGammaDensityObj::~InvGammaDensityObj() {
  // pass...
}

void
InvGammaDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
}

void
InvGammaDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform on arguments
  Tmin -= _shift; Tmax -= _shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::dinvgamma(t, _alpha, _beta);
  }
}

void
InvGammaDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform on arguments
  Tmin -= _shift; Tmax -= _shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::pinvgamma(t, _alpha, _beta);
  }
}

void
InvGammaDensityObj::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = RNG::invgamma(_alpha, _beta)+_shift;
  }
}

Density
InvGammaDensityObj::affine(double scale, double shift) const {
  return new GammaDensityObj(_alpha, scale*_beta, scale*_shift+shift);
}

void
InvGammaDensityObj::rangeEst(double alpha, double &a, double &b) const {
  a = stochbb::qinvgamma(alpha/2, _alpha, _beta);
  b = stochbb::qinvgamma(1-alpha/2, _alpha, _beta);
}

int
InvGammaDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // Compare uniform densities
  const InvGammaDensityObj *ogamma = dynamic_cast<const InvGammaDensityObj *>(&other);
  if (_alpha < ogamma->_alpha) { return -1; }
  else if (_alpha > ogamma->_alpha) { return 1; }
  if (_beta < ogamma->_beta) { return -1; }
  else if (_beta > ogamma->_beta) { return 1; }
  if (_shift < ogamma->_shift) { return -1; }
  else if (_shift > ogamma->_shift) { return 1; }
  return 0;
}

void
InvGammaDensityObj::print(std::ostream &stream) const {
  stream << "<InvGammaDensityObj alpha=" << _alpha << ", beta=" << _beta;
  if (_shift) { stream << ", shift=" << _shift; }
  stream << " #" << (void *)this << ">";
}


/* ********************************************************************************************* *
 * Implementation of WeibullDensityObj
 * ********************************************************************************************* */
WeibullDensityObj::WeibullDensityObj(double k, double lambda, double shift) throw (Error)
  : AtomicDensityObj(), _k(k), _lambda(lambda), _shift(shift)
{
  logDebug() << "Create WeibullDensity with k=" << _k << ", lambda=" << _lambda
             << ", shift=" << _shift << " #" << this << "." ;
  assume(_k > 0);
  assume(_lambda > 0);
}

WeibullDensityObj::~WeibullDensityObj() {
  // pass...
}

void
WeibullDensityObj::mark() {
  if (isMarked()) { return; }
  AtomicDensityObj::mark();
}

void
WeibullDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform on arguments
  Tmin -= _shift; Tmax -= _shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::dweibull(t, _k, _lambda);
  }
}

void
WeibullDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform on arguments
  Tmin -= _shift; Tmax -= _shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::pweibull(t, _k, _lambda);
  }
}

void
WeibullDensityObj::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = _lambda*std::pow(-std::log(1-RNG::unif()), 1./_k) + _shift;
  }
}

Density
WeibullDensityObj::affine(double scale, double shift) const {
  return new WeibullDensityObj(_k, scale*_lambda, scale*_shift+shift);
}

void
WeibullDensityObj::rangeEst(double alpha, double &a, double &b) const {
  a = stochbb::qweibull(alpha/2, _k, _lambda);
  b = stochbb::qweibull(1-alpha/2, _k, _lambda);
}

int
WeibullDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) {
    return res;
  }
  // Compare Weibull densities
  const WeibullDensityObj *owb = dynamic_cast<const WeibullDensityObj *>(&other);
  if (_k < owb->_k) { return -1; }
  else if (_k > owb->_k) { return 1; }
  if (_lambda < owb->_lambda) { return -1; }
  else if (_lambda > owb->_lambda) { return 1; }
  if (_shift < owb->_shift) { return -1; }
  else if (_shift > owb->_shift) { return 1; }
  return 0;
}

void
WeibullDensityObj::print(std::ostream &stream) const {
  stream << "<WeibullDensityObj k=" << _k << ", lambda=" << _lambda;
  if (_shift) { stream << ", shift=" << _shift; }
  stream << " #" << (void *)this << ">";
}


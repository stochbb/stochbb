#include "distribution.hh"
#include "math.hh"
#include "rng.hh"

using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of DistributionObj
 * ********************************************************************************************* */
DistributionObj::DistributionObj()
  : Object()
{
  // pass...
}

DistributionObj::~DistributionObj() {
  // pass...
}

void
DistributionObj::mark() {
  if (isMarked()) { return; }
  Object::mark();
}

/* ********************************************************************************************* *
 * Implementation of GenericAtomicDensityObj
 * ********************************************************************************************* */
GenericAtomicDensityObj::GenericAtomicDensityObj(DistributionObj *dist, Eigen::Ref<Eigen::VectorXd> params)
  : AtomicDensityObj(), _distribution(dist), _params(params)
{
  // pass...
}

GenericAtomicDensityObj::~GenericAtomicDensityObj() {
  // pass...
}

void
GenericAtomicDensityObj::mark() {
  if (isMarked()) { return; }
  AtomicDensityObj::mark();
  if (_distribution) { _distribution->mark(); }
}

void
GenericAtomicDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  double t = Tmin, dt=(Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out(i) = _distribution->pdf(t, _params);
  }
}

void
GenericAtomicDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  double t = Tmin, dt=(Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out(i) = _distribution->cdf(t, _params);
  }
}

void
GenericAtomicDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // Compare uniform densities
  const GenericAtomicDensityObj *ogen = dynamic_cast<const GenericAtomicDensityObj *>(&other);
  // Compare by distribution type
  if (_a < ounif->_a) { return -1; }
  else if (_a > ounif->_a) { return 1; }
  if (_b < ounif->_b) { return -1; }
  else if (_b > ounif->_b) { return 1; }
  return 0;

}


/* ********************************************************************************************* *
 * Implementation of UniformDistObj
 * ********************************************************************************************* */
UniformDistObj::UniformDistObj()
  : DistributionObj()
{
  // pass...
}

UniformDistObj::~UniformDistObj() {
  // pass...
}

void
UniformDistObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
UniformDistObj::nParams() const {
  return 2;
}

double
UniformDistObj::pdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  if ((x < a) || (x > b)) { return 0; }
  return 1./(b-a);
}

double
UniformDistObj::cdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  if (x < a) { return 0.0; }
  else if (x > b) { return 1.0; }
  return (x-a)/(b-a);
}

double
UniformDistObj::quantile(double p, const Eigen::Ref<Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  double d = (b-a)*p/2;
  return a+d;
}

double
UniformDistObj::sample(const Eigen::Ref<Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  return RNG::unif()*(b-a) + a;
}


/* ********************************************************************************************* *
 * Implementation of NormalDistObj
 * ********************************************************************************************* */
NormalDistObj::NormalDistObj()
  : DistributionObj()
{
  // pass...
}

NormalDistObj::~NormalDistObj() {
  // pass...
}

void
NormalDistObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
NormalDistObj::nParams() const {
  return 2;
}

double
NormalDistObj::pdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  return stochbb::dnorm((x-mu)/sigma);
}

double
NormalDistObj::cdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  return stochbb::pnorm((x-mu)/sigma);
}

double
NormalDistObj::quantile(double p, const Eigen::Ref<Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  return stochbb::qnorm(p)*sigma + mu;
}

double
NormalDistObj::sample(const Eigen::Ref<Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  return RNG::norm()*sigma + mu;
}


/* ********************************************************************************************* *
 * Implementation of GammaDistObj
 * ********************************************************************************************* */
GammaDistObj::GammaDistObj()
  : DistributionObj()
{
  // pass...
}

GammaDistObj::~GammaDistObj() {
  // pass...
}

void
GammaDistObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
GammaDistObj::nParams() const {
  return 3;
}

double
GammaDistObj::pdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  return stochbb::dgamma(x-shift, k, theta);
}

double
GammaDistObj::cdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  return stochbb::pgamma(x-shift, k, theta);
}

double
GammaDistObj::quantile(double p, const Eigen::Ref<Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  return stochbb::qgamma(p, k, theta) + shift;
}

double
GammaDistObj::sample(const Eigen::Ref<Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  return RNG::gamma(k, theta) + shift;
}


/* ********************************************************************************************* *
 * Implementation of InvGammaDistObj
 * ********************************************************************************************* */
InvGammaDistObj::InvGammaDistObj()
  : DistributionObj()
{
  // pass...
}

InvGammaDistObj::~InvGammaDistObj() {
  // pass...
}

void
InvGammaDistObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
InvGammaDistObj::nParams() const {
  return 3;
}

double
InvGammaDistObj::pdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  return stochbb::dinvgamma(x-shift, alpha, beta);
}

double
InvGammaDistObj::cdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  return stochbb::pinvgamma(x-shift, alpha, beta);
}

double
InvGammaDistObj::quantile(double p, const Eigen::Ref<Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  return stochbb::qinvgamma(p, alpha, beta) + shift;
}

double
InvGammaDistObj::sample(const Eigen::Ref<Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  return RNG::invgamma(alpha, beta) + shift;
}


/* ********************************************************************************************* *
 * Implementation of WeibullDistObj
 * ********************************************************************************************* */
WeibullDistObj::WeibullDistObj()
  : DistributionObj()
{
  // pass...
}

WeibullDistObj::~WeibullDistObj() {
  // pass...
}

void
WeibullDistObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
WeibullDistObj::nParams() const {
  return 3;
}

double
WeibullDistObj::pdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double k=params[0], lambda=params[1], shift=params[2];
  return stochbb::dweibull(x-shift, k, lambda);
}

double
WeibullDistObj::cdf(double x, const Eigen::Ref<Eigen::VectorXd> params) const {
  double k=params[0], lambda=params[1], shift=params[2];
  return stochbb::pweibull(x-shift, k, lambda);
}

double
WeibullDistObj::quantile(double p, const Eigen::Ref<Eigen::VectorXd> params) const {
  double k=params[0], lambda=params[1], shift=params[2];
  return stochbb::qweibull(p, k, lambda) + shift;
}

double
WeibullDistObj::sample(const Eigen::Ref<Eigen::VectorXd> params) const {
  double k=params[0], lambda=params[1], shift=params[2];
  return RNG::invgamma(k, lambda) + shift;
}

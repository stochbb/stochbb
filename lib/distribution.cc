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

int
DistributionObj::compare(const DistributionObj &other) const {
  // Same object -> equal
  if (this == &other) { return 0; }
  // compare by type
  if (typeid(*this).before(typeid(other))) { return -1; }
  else if (typeid(other).before(typeid(*this))) { return -1; }
  return 0;
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

Density
GenericAtomicDensityObj::affine(double scale, double shift) const {
  // copy parameters
  Eigen::VectorXd params = _params;
  _distribution->affine(scale, shift, params);
  return new GenericAtomicDensityObj(_distribution, params);
}

void
GenericAtomicDensityObj::rangeEst(double alpha, double &a, double &b) const {

}

void
GenericAtomicDensityObj::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  for (int i=0; i<out.size(); i++) {
    out[i] = _distribution->sample(_params);
  }
}

int
GenericAtomicDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // Compare atomic densities
  const GenericAtomicDensityObj *ogen = dynamic_cast<const GenericAtomicDensityObj *>(&other);
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
GenericAtomicDensityObj::print(std::ostream &stream) const {
  stream << "<AtomicDensity ";
  _distribution->print(stream);
  stream << " @" << _params.transpose() << " #" << this << ">";
}


/* ********************************************************************************************* *
 * Implementation of UniformDistObj
 * ********************************************************************************************* */
UniformDistributionObj::UniformDistributionObj()
  : DistributionObj()
{
  // pass...
}

UniformDistributionObj::~UniformDistributionObj() {
  // pass...
}

void
UniformDistributionObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
UniformDistributionObj::nParams() const {
  return 2;
}

double
UniformDistributionObj::pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  if ((x < a) || (x > b)) { return 0; }
  return 1./(b-a);
}

double
UniformDistributionObj::cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  if (x < a) { return 0.0; }
  else if (x > b) { return 1.0; }
  return (x-a)/(b-a);
}

void
UniformDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  double d = (b-a)*p/2;
  lower = a+d; upper = b-d;
}

void
UniformDistributionObj::affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const {
  double a = params[0], b = params[1];
  a = scale*a + shift;
  b = scale*b + shift;
  params[0] = std::min(a,b);
  params[1] = std::max(a,b);
}

void
UniformDistributionObj::affine(double scale, double shift, std::vector<DensityObj *> &params) const {
  if (scale>0) {
    params[0] = *(params[0]->affine(scale, shift));
    params[1] = *(params[1]->affine(scale, shift));
  } else {
    params[0] = *(params[0]->affine(scale, shift));
    params[1] = *(params[1]->affine(scale, shift));
    std::swap(params[0], params[1]);
  }
}

double
UniformDistributionObj::sample(const Eigen::Ref<const Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  return RNG::unif()*(b-a) + a;
}


/* ********************************************************************************************* *
 * Implementation of NormalDistObj
 * ********************************************************************************************* */
NormalDistributionObj::NormalDistributionObj()
  : DistributionObj()
{
  // pass...
}

NormalDistributionObj::~NormalDistributionObj() {
  // pass...
}

void
NormalDistributionObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
NormalDistributionObj::nParams() const {
  return 2;
}

double
NormalDistributionObj::pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  return stochbb::dnorm((x-mu)/sigma);
}

double
NormalDistributionObj::cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  return stochbb::pnorm((x-mu)/sigma);
}

void NormalDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  lower=mu-stochbb::qnorm(p)*sigma; upper=mu+stochbb::qnorm(p)*sigma;
}

double
NormalDistributionObj::sample(const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  return RNG::norm()*sigma + mu;
}


/* ********************************************************************************************* *
 * Implementation of GammaDistObj
 * ********************************************************************************************* */
GammaDistributionObj::GammaDistributionObj()
  : DistributionObj()
{
  // pass...
}

GammaDistributionObj::~GammaDistributionObj() {
  // pass...
}

void
GammaDistributionObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
GammaDistributionObj::nParams() const {
  return 3;
}

double
GammaDistributionObj::pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  return stochbb::dgamma(x-shift, k, theta);
}

double
GammaDistributionObj::cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  return stochbb::pgamma(x-shift, k, theta);
}

void GammaDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  lower=shift; upper = stochbb::qgamma(p, k, theta) + shift;
}

double
GammaDistributionObj::sample(const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  return RNG::gamma(k, theta) + shift;
}


/* ********************************************************************************************* *
 * Implementation of InvGammaDistObj
 * ********************************************************************************************* */
InvGammaDistributionObj::InvGammaDistributionObj()
  : DistributionObj()
{
  // pass...
}

InvGammaDistributionObj::~InvGammaDistributionObj() {
  // pass...
}

void
InvGammaDistributionObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
InvGammaDistributionObj::nParams() const {
  return 3;
}

double
InvGammaDistributionObj::pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  return stochbb::dinvgamma(x-shift, alpha, beta);
}

double
InvGammaDistributionObj::cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  return stochbb::pinvgamma(x-shift, alpha, beta);
}

void InvGammaDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  lower=shift; upper = stochbb::qinvgamma(p, alpha, beta) + shift;
}

double
InvGammaDistributionObj::sample(const Eigen::Ref<const Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  return RNG::invgamma(alpha, beta) + shift;
}


/* ********************************************************************************************* *
 * Implementation of WeibullDistObj
 * ********************************************************************************************* */
WeibullDistributionObj::WeibullDistributionObj()
  : DistributionObj()
{
  // pass...
}

WeibullDistributionObj::~WeibullDistributionObj() {
  // pass...
}

void
WeibullDistributionObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
WeibullDistributionObj::nParams() const {
  return 3;
}

double
WeibullDistributionObj::pdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], lambda=params[1], shift=params[2];
  return stochbb::dweibull(x-shift, k, lambda);
}

double
WeibullDistributionObj::cdf(double x, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], lambda=params[1], shift=params[2];
  return stochbb::pweibull(x-shift, k, lambda);
}

void WeibullDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], lambda=params[1], shift=params[2];
  lower = shift; upper = stochbb::qweibull(p, k, lambda) + shift;
}

double
WeibullDistributionObj::sample(const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], lambda=params[1], shift=params[2];
  return RNG::invgamma(k, lambda) + shift;
}


/* ********************************************************************************************* *
 * Implementation of GenericCompoundDensityObj
 * ********************************************************************************************* */
GenericCompoundDensityObj::GenericCompoundDensityObj(DistributionObj *dist, const std::vector<DensityObj *> &params)
  : DensityObj(), _distribution(dist), _parameters(params)
{
  // pass...
}

GenericCompoundDensityObj::~GenericCompoundDensityObj() {
  // pass...
}

void
GenericCompoundDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  _distribution->mark();
  for (size_t i=0; i<_parameters.size(); i++) {
    _parameters[i]->mark();
  }
}

void
GenericCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  size_t Nstep = 100;
  double df = 1;

  // Get parameter values and PDFs
  Eigen::MatrixXd PDFs(Nstep, _parameters.size());
  Eigen::MatrixXd params(Nstep, _parameters.size());
  for (size_t i=0; i<_parameters.size(); i++) {
    double a, b;
    _parameters[i]->rangeEst(0.01, a, b);

    double x = a, dx = (b-a)/Nstep;
    df *= (b-a)/Nstep;

    for (size_t j=0; j<Nstep; j++, x+=dx)
      params(j,i) = x;

    _parameters[i]->eval(a,b, PDFs.col(i));
  }

  out.setZero();
  size_t N = std::pow(Nstep, _parameters.size());
  Eigen::VectorXd param(_parameters.size());
  std::vector<size_t> idxs(0, _parameters.size());
  double x=Tmin, dx=(Tmax-Tmin)/out.size();

  // For each value in [Tmin, Tmax):
  for (int i=0; i<out.size(); i++, x+=dx) {
    // "Integrate" over parameter space
    for (size_t j=0; j<N; j++) {
      // Get parameter indices
      _to_param_indices(j, Nstep, idxs);
      // Get parameter vector and prod of parameter PDFs
      double dp = 1;
      for (size_t k=0; k<_parameters.size(); k++) {
        param[k] = params(idxs[k], k);
        dp *= PDFs(idxs[k], k);
      }
      out(i) += _distribution->pdf(x, param)*dp*df;
    }
  }
}

void
GenericCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  size_t Nstep = 100;
  double df = 1;

  // Get parameter values and PDFs
  Eigen::MatrixXd PDFs(Nstep, _parameters.size());
  Eigen::MatrixXd params(Nstep, _parameters.size());
  for (size_t i=0; i<_parameters.size(); i++) {
    double a, b;
    _parameters[i]->rangeEst(0.01, a, b);

    double x = a, dx = (b-a)/Nstep;
    df *= (b-a)/Nstep;

    for (size_t j=0; j<Nstep; j++, x+=dx)
      params(j,i) = x;

    _parameters[i]->eval(a,b, PDFs.col(i));
  }

  out.setZero();
  size_t N = std::pow(Nstep, _parameters.size());
  Eigen::VectorXd param(_parameters.size());
  std::vector<size_t> idxs(0, _parameters.size());
  double x=Tmin, dx=(Tmax-Tmin)/out.size();

  // For each value in [Tmin, Tmax):
  for (int i=0; i<out.size(); i++, x+=dx) {
    // "Integrate" over parameter space
    for (size_t j=0; j<N; j++) {
      // Get parameter indices
      _to_param_indices(j, Nstep, idxs);
      // Get parameter vector and prod of parameter PDFs
      double dp = 1;
      for (size_t k=0; k<_parameters.size(); k++) {
        param[k] = params(idxs[k], k);
        dp *= PDFs(idxs[k], k);
      }
      out(i) += _distribution->cdf(x, param)*dp*df;
    }
  }
}

void
GenericCompoundDensityObj::_to_param_indices(size_t i, size_t N, std::vector<size_t> &idxs) const {
  for (size_t j=0; j<idxs.size(); j++) {
    idxs[j] = i%N; i /= N;
  }
}

Density
GenericCompoundDensityObj::affine(double scale, double shift) const {
  GenericCompoundDensityObj *res = new GenericCompoundDensityObj(_distribution, _parameters);
  _distribution->affine(scale, shift, res->_parameters);
  return res;
}

void
GenericCompoundDensityObj::rangeEst(double alpha, double &a, double &b) const {
  Eigen::VectorXd A(_parameters.size()), B(_parameters.size());
  for (size_t i=0; i<_parameters.size(); i++) {
    _parameters[i]->rangeEst(alpha, A[i], B[i]);
  }
  double tmp;
  _distribution->quantile(a, tmp, alpha, A);
  _distribution->quantile(tmp, b, alpha, B);
}

int
GenericCompoundDensityObj::compare(const DensityObj &other) const {
  // Compare by density type
  if (int res = DensityObj::compare(other))
    return res;
  // same density class -> cast
  const GenericCompoundDensityObj &cother = static_cast<const GenericCompoundDensityObj &>(other);
  // Compare by distribution
  if (int res = _distribution->compare(*cother._distribution))
    return res;
  // Compare by parameter densities
  for (size_t i=0; i<_parameters.size(); i++) {
    if (int res = _parameters[i]->compare(*cother._parameters[i]))
      return res;
  }
  return 0;
}

void
GenericCompoundDensityObj::print(std::ostream &stream) const {
  stream << "<CompoundDensity of "; _distribution->print(stream);
  stream << " with [";
  _parameters[0]->print(stream);
  for (size_t i=1; i<_parameters.size(); i++) {
    stream << ", ";
    _parameters[i]->print(stream);
  }
  stream << "] #" << this << ">";
}

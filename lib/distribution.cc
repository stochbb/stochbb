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

GenericAtomicDensityObj::GenericAtomicDensityObj(const Distribution &dist, Eigen::Ref<Eigen::VectorXd> params)
  : AtomicDensityObj(), _distribution(*dist), _params(params)
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

Distribution
GenericAtomicDensityObj::distribution() const {
  _distribution->ref();
  return _distribution;
}

size_t
GenericAtomicDensityObj::nParams() const {
  return _distribution->nParams();
}

double
GenericAtomicDensityObj::parameter(size_t i) const {
  return _params[i];
}

void
GenericAtomicDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  _distribution->pdf(Tmin, Tmax, out, _params);
}

void
GenericAtomicDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  _distribution->cdf(Tmin, Tmax, out, _params);
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
  /// @bug Implement!
}

void
GenericAtomicDensityObj::sample(Eigen::Ref<Eigen::VectorXd> out) const {
  _distribution->sample(out, _params);
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
 * Implementation of DeltaDistObj
 * ********************************************************************************************* */
DeltaDistributionObj::DeltaDistributionObj()
  : DistributionObj()
{
  // pass...
}

DeltaDistributionObj::~DeltaDistributionObj() {
  // pass...
}

void
DeltaDistributionObj::mark() {
  if (isMarked()) { return; }
  DistributionObj::mark();
}

size_t
DeltaDistributionObj::nParams() const {
  return 1;
}

void
DeltaDistributionObj::pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
                          const Eigen::Ref<const Eigen::VectorXd> params) const
{
  out.setZero();
  if ((params[0]<Tmin) || (params[0]>Tmax)) { return; }
  double dt = (Tmax-Tmin)/out.size();
  size_t idx = (params[0]-Tmin)/dt;
  out[idx] = 1./dt;
}

void
DeltaDistributionObj::cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
                          const Eigen::Ref<const Eigen::VectorXd> params) const
{
  out.setZero();
  if ((params[0]<Tmin) || (params[0]>Tmax)) { return; }
  double dt = (Tmax-Tmin)/out.size();
  size_t idx = (params[0]-Tmin)/dt;
  out.tail(out.size()-idx).setConstant(1);
}

void
DeltaDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  lower = upper = params[0];
}

void
DeltaDistributionObj::affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const {
  params[0] = scale*params[0] + shift;
}

void
DeltaDistributionObj::affine(double scale, double shift, std::vector<DensityObj *> &params) const {
  params[0] = *(params[0]->affine(scale, shift));
}

void
DeltaDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  out.setConstant(params[0]);
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

void
UniformDistributionObj::pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = ((t >= a) && (t <= b)) ? 1./(b-a) : 0.0;
  }
}

void
UniformDistributionObj::cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    if (t<a) { out[i] = 0.0; }
    else if (t>b) { out[i] = 1.0; }
    else { out[i] = (t-a)/(b-a); }
  }
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

void
UniformDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  for (int i=0; i<out.size(); i++) {
   out(i) = RNG::unif()*(b-a) + a;
  }
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

void NormalDistributionObj::pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::dnorm((t-mu)/sigma)/sigma;
  }
}

void NormalDistributionObj::cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::pnorm((t-mu)/sigma);
  }
}

void NormalDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  lower=mu-stochbb::qnorm(p)*sigma; upper=mu+stochbb::qnorm(p)*sigma;
}

void NormalDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  for (int i=0; i<out.size(); i++) {
    out(i) = RNG::norm()*sigma + mu;
  }
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

void GammaDistributionObj::pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  Tmin -= shift; Tmax -= shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::dgamma(t, k, theta);
  }
}

void GammaDistributionObj::cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  Tmin -= shift; Tmax -= shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::pgamma(t, k, theta);
  }
}

void GammaDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  lower=shift; upper = stochbb::qgamma(p, k, theta) + shift;
}

void
GammaDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  for (int i=0; i<out.size(); i++) {
    out(i) = RNG::gamma(k, theta) + shift;
  }
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

void InvGammaDistributionObj::pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  Tmin -= shift; Tmax -= shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::dinvgamma(t, alpha, beta);
  }
}

void InvGammaDistributionObj::cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  Tmin -= shift; Tmax -= shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::pinvgamma(t, alpha, beta);
  }
}

void InvGammaDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  lower=shift; upper = stochbb::qinvgamma(p, alpha, beta) + shift;
}

void InvGammaDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double alpha=params[0], beta=params[1], shift=params[2];
  for (int i=0; i<out.size(); i++) {
    out(i) = RNG::invgamma(alpha, beta) + shift;
  }
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

void
WeibullDistributionObj::pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
                            const Eigen::Ref<const Eigen::VectorXd> params) const
{
  double k=params[0], lambda=params[1], shift=params[2];
  // Apply affine transform on arguments
  Tmin -= shift; Tmax -= shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::dweibull(t, k, lambda);
  }
}

void
WeibullDistributionObj::cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
                            const Eigen::Ref<const Eigen::VectorXd> params) const
{
  double k=params[0], lambda=params[1], shift=params[2];
  Tmin -= shift; Tmax -= shift;
  double t = Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out[i] = stochbb::pweibull(t, k, lambda);
  }
}

void
WeibullDistributionObj::quantile(double &lower, double &upper, double p,
                                 const Eigen::Ref<const Eigen::VectorXd> params) const
{
  double k=params[0], lambda=params[1], shift=params[2];
  lower = shift; upper = stochbb::qweibull(p, k, lambda) + shift;
}

void
WeibullDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], lambda=params[1], shift=params[2];
  for (int i=0; i<out.size(); i++) {
    out(i) = RNG::invgamma(k, lambda) + shift;
  }
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

Distribution
GenericCompoundDensityObj::distribution() const {
  _distribution->ref();
  return _distribution;
}

size_t
GenericCompoundDensityObj::nParams() const {
  return _distribution->nParams();
}

Density
GenericCompoundDensityObj::parameter(size_t i) const {
  _parameters[i]->ref();
  return _parameters[i];
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
  double dx=(Tmax-Tmin)/out.size();

  // For each value in [Tmin, Tmax):
#pragma omp for
  for (int i=0; i<out.size(); i++) {
    double x = Tmin + i*dx;
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
      _distribution->pdf(x, x+dx, out.segment<1>(i), param);
      out(i) *= dp*df;
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
  double dx=(Tmax-Tmin)/out.size();

  // For each value in [Tmin, Tmax):
#pragma omp for
  for (int i=0; i<out.size(); i++) {
    double x = Tmin + i*dx;
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
      _distribution->cdf(x, x+dx, out.segment<1>(i), param);
      out(i) *= dp*df;
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

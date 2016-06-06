#include "distribution.hh"
#include "math.hh"
#include "rng.hh"
#include <iostream>
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
  if ((params[0]<Tmin) || (params[0]>=Tmax)) { return; }
  double dt = (Tmax-Tmin)/out.size();
  int idx = (params[0]-Tmin)/dt;
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
DeltaDistributionObj::quantile(double &lower, double &upper, double p,
                               const Eigen::Ref<const Eigen::VectorXd> params) const {
  lower = params[0]-1e-6; upper = params[0]+1e-6;
}

void
DeltaDistributionObj::affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const {
  params[0] = scale*params[0] + shift;
}

void
DeltaDistributionObj::affine(double scale, double shift, std::vector<Density> &params) const {
  params[0] = params[0].affine(scale, shift);
}

void
DeltaDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out,
                             const Eigen::Ref<const Eigen::VectorXd> params) const {
  out.setConstant(params[0]);
}

void
DeltaDistributionObj::print(std::ostream &stream) const {
  stream << "<DeltaDistribution #" << this << ">";
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
UniformDistributionObj::affine(double scale, double shift, std::vector<Density> &params) const {
  if (scale>0) {
    params[0] = params[0].affine(scale, shift);
    params[1] = params[1].affine(scale, shift);
  } else {
    params[0] = params[0].affine(scale, shift);
    params[1] = params[1].affine(scale, shift);
    std::swap(params[0], params[1]);
  }
}

void
UniformDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double a=params[0], b=params[1];
  RNG &rng = RNG::get();
  std::uniform_real_distribution<double> sampler(a, b);
  for (int i=0; i<out.size(); i++) {
   out(i) = sampler(rng);
  }
}

void
UniformDistributionObj::print(std::ostream &stream) const {
  stream << "<UniformDistribution #" << this << ">";
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

void
NormalDistributionObj::affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const {
  double mu = params[0], sigma = params[1];
  params[0] = scale*mu + shift;
  params[1] = scale*sigma;
}

void
NormalDistributionObj::affine(double scale, double shift, std::vector<Density> &params) const {
  Density mu = params[0], sigma = params[1];
  params[0] = mu.affine(scale, shift);
  params[1] = sigma.affine(scale,0);
}

void NormalDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  lower=mu+stochbb::qnorm(p)*sigma; upper=mu+stochbb::qnorm(1-p)*sigma;
}

void NormalDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double mu=params[0], sigma=params[1];
  RNG &rng = RNG::get();
  std::normal_distribution<double> sampler(mu, sigma);
  for (int i=0; i<out.size(); i++) {
    out(i) = sampler(rng);
  }
}

void
NormalDistributionObj::print(std::ostream &stream) const {
  stream << "<NormalDistribution #" << this << ">";
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

void
GammaDistributionObj::affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const {
  double theta = params[1], s = params[2];
  params[1] = scale*theta;
  params[2] = scale*s+shift;
}

void
GammaDistributionObj::affine(double scale, double shift, std::vector<Density> &params) const {
  Density theta = params[1], s = params[2];
  params[1] = theta.affine(scale, 0);
  params[2] = s.affine(scale,shift);
}

void GammaDistributionObj::quantile(double &lower, double &upper, double p, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  lower=shift; upper = stochbb::qgamma(p, k, theta) + shift;
}

void
GammaDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd> params) const {
  double k=params[0], theta=params[1], shift=params[2];
  RNG &rng = RNG::get();
  std::gamma_distribution<double> sampler(k, theta);
  for (int i=0; i<out.size(); i++) {
    out(i) = sampler(rng) + shift;
  }
}

void
GammaDistributionObj::print(std::ostream &stream) const {
  stream << "<GammaDistribution #" << this << ">";
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

void
InvGammaDistributionObj::affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const {
  double beta = params[1], s = params[2];
  params[1] = scale*beta;
  params[2] = scale*s+shift;
}

void
InvGammaDistributionObj::affine(double scale, double shift, std::vector<Density> &params) const {
  Density beta = params[1], s = params[2];
  params[1] = beta.affine(scale, 0);
  params[2] = s.affine(scale,shift);
}

void
InvGammaDistributionObj::quantile(double &lower, double &upper, double p,
                                  const Eigen::Ref<const Eigen::VectorXd> params) const
{
  double alpha=params[0], beta=params[1], shift=params[2];
  lower=shift; upper = stochbb::qinvgamma(p, alpha, beta) + shift;
}

void
InvGammaDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out,
                                const Eigen::Ref<const Eigen::VectorXd> params) const
{
  double alpha=params[0], beta=params[1], shift=params[2];
  RNG &rng=RNG::get();
  std::inverse_gamma_distribution<double> sampler(alpha, 1./beta);
  for (int i=0; i<out.size(); i++) {
    out(i) =  sampler(rng) + shift;
  }
}

void
InvGammaDistributionObj::print(std::ostream &stream) const {
  stream << "<InvGammaDistribution #" << this << ">";
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
WeibullDistributionObj::affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const {
  double lambda = params[1], s = params[2];
  params[1] = scale*lambda;
  params[2] = scale*s+shift;
}

void
WeibullDistributionObj::affine(double scale, double shift, std::vector<Density> &params) const {
  Density lambda = params[1], s = params[2];
  params[1] = lambda.affine(scale, 0);
  params[2] = s.affine(scale,shift);
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
  RNG &rng = RNG::get();
  std::weibull_distribution<double> sampler(k, lambda);
  for (int i=0; i<out.size(); i++) {
    out(i) = sampler(rng) + shift;
  }
}

void
WeibullDistributionObj::print(std::ostream &stream) const {
  stream << "<WeibullDistribution #" << this << ">";
}


/* ********************************************************************************************* *
 * Implementation of StudtDistObj
 * ********************************************************************************************* */
StudtDistributionObj::StudtDistributionObj()
  : DistributionObj()
{
  // pass...
}

StudtDistributionObj::~StudtDistributionObj() {
  // pass...
}

size_t
StudtDistributionObj::nParams() const {
  return 3;
}

void
StudtDistributionObj::pdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
                          const Eigen::Ref<const Eigen::VectorXd> params) const
{
  double nu=params[0], sigma=params[1], mu=params[2];
  Tmin = (Tmin - mu)/sigma;
  Tmax = (Tmax - mu)/sigma;
  double t=Tmin, dt=(Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    out(i) = std::exp(
          std::lgamma((nu+1)/2)-std::lgamma(nu/2) - 0.5*std::log(M_PI*nu) - std::log(sigma)
          -(nu+1)*std::log(1+t*t/nu)/2);
  }
}

void
StudtDistributionObj::cdf(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out,
                          const Eigen::Ref<const Eigen::VectorXd> params) const
{
  double nu=params[0], sigma=params[1], mu=params[2];
  Tmin = (Tmin - mu)/sigma;
  Tmax = (Tmax - mu)/sigma;
  double t=Tmin, dt=(Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    double x = nu/(t*t+nu);
    if (0 >= t) {
      out(i) = 0.5*stochbb::betai(nu/2, 0.5, x);
    } else {
      out(i) = 1-0.5*stochbb::betai(nu/2, 0.5, x);
    }
  }
}

void
StudtDistributionObj::quantile(double &lower, double &upper, double p,
                               const Eigen::Ref<const Eigen::VectorXd> params) const
{
  double nu=params[0], sigma=params[1], mu=params[2];
  // get inverse reg. incomplete beta function
  upper = invbetai(nu/2, 0.5, 1-p/2);
  // map
  upper = std::sqrt(nu/upper - nu);
  // symm.
  lower = -upper;
  // scale/shift
  upper = sigma*upper+mu;
  lower = sigma*lower+mu;
}

void
StudtDistributionObj::affine(double scale, double shift, Eigen::Ref<Eigen::VectorXd> params) const {
  params[1] *= scale;
  params[2] = scale*params[2] + shift;
}

void
StudtDistributionObj::affine(double scale, double shift, std::vector<Density> &params) const {
  Density sigma = params[1], mu = params[2];
  params[1] = sigma.affine(scale, 0);
  params[2] = mu.affine(scale,shift);
}

void
StudtDistributionObj::sample(Eigen::Ref<Eigen::VectorXd> out,
                             const Eigen::Ref<const Eigen::VectorXd> params) const
{
  double df=params[0], sigma=params[1], mu = params[2];
  RNG &rng = RNG::get();
  std::student_t_distribution<double> sampler(df);
  for (int i=0; i<out.size(); i++) {
    out(i) = sigma*sampler(rng) + mu;
  }
}

void
StudtDistributionObj::print(std::ostream &stream) const {
  stream << "<Student's t-Distribution #" << this << ">";
}

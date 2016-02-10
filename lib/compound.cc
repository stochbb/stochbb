#include "compound.hh"
#include "operators.hh"
#include "math.hh"

using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of CompoundObj
 * ********************************************************************************************* */
CompoundObj::CompoundObj(const std::vector<Var> &vars, const std::string &name)
  : DerivedVarObj(vars, name)
{
  // pass...
}

void
CompoundObj::mark() {
  if (isMarked()) { return; }
  DerivedVarObj::mark();
}

Compound
CompoundObj::norm(const Var &mu, const Var &sigma, const std::string &name) {
  return new NormalCompoundObj(mu, sigma, name);
}

Compound
CompoundObj::gamma(const Var &k, const Var &theta, const std::string &name) {
  return new GammaCompoundObj(k, theta, name);
}


/* ********************************************************************************************* *
 * Implementation of NormalCompoundDensityObj
 * ********************************************************************************************* */
NormalCompoundDensityObj::NormalCompoundDensityObj(DensityObj *mu, DensityObj *sigma)
  : DensityObj(), _mu(mu), _sigma(sigma), _muMin(0), _ddMu(0), _dmu(),
    _sigMin(0), _ddSig(0), _dsigma()
{
  // Prepare integration
  _init_int();
}

NormalCompoundDensityObj::NormalCompoundDensityObj(const Var &mu, const Var &sigma) throw (AssumptionError)
  : DensityObj(), _mu(*mu->density()), _sigma(*sigma->density()), _muMin(0), _ddMu(0), _dmu(),
    _sigMin(0), _ddSig(0), _dsigma()
{
  // Test for independence
  if (! mu.mutuallyIndep(sigma)) {
    AssumptionError err;
    err << "Cannot construct normal-compound density: parameters are not independent.";
    throw err;
  }

  // Prepare integration
  _init_int();
}

void
NormalCompoundDensityObj::_init_int() {
  double muMax, sigMax;
  // Get intergration range for mu
  _mu->rangeEst(0.01, _muMin, muMax);
  // Get intergration range for sigma
  _sigma->rangeEst(0.01, _sigMin, sigMax);

  // Check if mu is a delta distribution
  if (DeltaDensityObj *delta_mu = dynamic_cast<DeltaDensityObj *>(_mu)) {
    _dmu.resize(1); _dmu[0] = 1; _ddMu = 1; _muMin = delta_mu->delay();
  } else {
    _dmu.resize(100);
    _ddMu = (muMax-_muMin)/100;
    _mu->eval(_muMin, muMax, _dmu);
  }

  // Check if sigma is a delta distribution
  if (DeltaDensityObj *delta_sig = dynamic_cast<DeltaDensityObj *>(_sigma)) {
    _dsigma.resize(1); _dsigma[0] = 1; _ddSig = 1; _sigMin = delta_sig->delay();
  } else {
    _dsigma.resize(100);
    _ddSig = (sigMax-_sigMin)/100;
    _sigma->eval(_sigMin, sigMax, _dsigma);
  }
}

NormalCompoundDensityObj::~NormalCompoundDensityObj() {
  // pass...
}

void
NormalCompoundDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  _mu->mark();
  _sigma->mark();
}

void
NormalCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double mu = _muMin;
    for (int k=0; k<_dmu.size(); k++, mu+=_ddMu) {
      double sigma = _sigMin;
      for (int l=0; l<_dsigma.size(); l++, sigma+=_ddSig) {
        out(i) += _ddMu*_ddSig * _dmu[k] * _dsigma[l] *
            std::exp( -(t-mu)*(t-mu)/(2*sigma*sigma)) / (std::sqrt(2*M_PI)*sigma);
      }
    }
  }
}


void
NormalCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double mu = _muMin;
    for (int k=0; k<_dmu.size(); k++, mu+=_ddMu) {
      double sigma = _sigMin;
      for (int l=0; l<_dsigma.size(); l++, sigma+=_ddSig) {
        out(i) += _ddMu*_dmu[k] * _ddSig*_dsigma[l] *
            0.5*(1+std::erf((t-mu)/(sigma*std::sqrt(2))));
      }
    }
  }
}

Density
NormalCompoundDensityObj::affine(double scale, double shift) const {
  return new NormalCompoundDensityObj(*_mu->affine(scale, shift), *_sigma->affine(scale, 0));
}

void
NormalCompoundDensityObj::rangeEst(double alpha, double &a, double &b) const {
  // Get quantiles from parameter distr.
  double a_mu, b_mu, a_sig, b_sig;
  _mu->rangeEst(alpha, a_mu, b_mu);
  _sigma->rangeEst(alpha, a_sig, b_sig);

  // Get max quantiles
  a = a_mu+qnorm(alpha)*b_sig;
  b = b_mu-qnorm(1-alpha)*b_sig;
}



/* ********************************************************************************************* *
 * Implementation of NormalCompoundObj
 * ********************************************************************************************* */
NormalCompoundObj::NormalCompoundObj(const Var &mu, const Var &sigma, const std::string &name)
  : CompoundObj(std::vector<Var> {mu, sigma}, name), _density(0)
{
  // Test for independence
  if (! mu.mutuallyIndep(sigma)) {
    AssumptionError err;
    err << "Cannot construct normal-compound parameters are not independent.";
    throw err;
  }

  _density = new NormalCompoundDensityObj(mu, sigma);
  _density->unref();
}

void
NormalCompoundObj::mark() {
  if (isMarked()) { return; }
  CompoundObj::mark();
  if (_density) { _density->mark(); }
}

Density
NormalCompoundObj::density() {
  _density->mark(); return _density;
}


/* ********************************************************************************************* *
 * Implementation of GammaCompoundDensityObj
 * ********************************************************************************************* */
GammaCompoundDensityObj::GammaCompoundDensityObj(DensityObj *k, DensityObj *theta, double shift)
  : DensityObj(), _k(k), _theta(theta), _shift(shift), _kMin(0), _ddK(0), _dk(),
    _thetaMin(0), _ddTheta(0), _dtheta()
{
  _init_int();
}

GammaCompoundDensityObj::GammaCompoundDensityObj(const Var &k, const Var &theta, double shift) throw (AssumptionError)
  : DensityObj(), _k(*k->density()), _theta(*theta->density()), _shift(shift), _kMin(0), _ddK(0), _dk(),
    _thetaMin(0), _ddTheta(0), _dtheta()
{
  // Test for independence
  if (! k.mutuallyIndep(theta)) {
    AssumptionError err;
    err << "Cannot construct gamma-compound parameters are not independent.";
    throw err;
  }

  _init_int();
}

void
GammaCompoundDensityObj::_init_int() {
  double kMax, thetaMax;
  // Get intergration range for k
  _k->rangeEst(0.01, _kMin, kMax);
  // Get intergration range for theta
  _theta->rangeEst(0.01, _thetaMin, thetaMax);

  // Check if theta is a delta distribution
  if (DeltaDensityObj *delta_k = dynamic_cast<DeltaDensityObj *>(_k)) {
    _dk.resize(1); _dk[0] = 1; _ddK = 1; _kMin = delta_k->delay();
  } else {
    _dk.resize(100);
    _ddK = (kMax-_kMin)/100;
    _k->eval(_kMin, kMax, _dk);
  }

  // Check if sigma is a delta distribution
  if (DeltaDensityObj *delta_theta = dynamic_cast<DeltaDensityObj *>(_theta)) {
    _dtheta.resize(1); _dtheta[0] = 1; _ddTheta = 1; _thetaMin = delta_theta->delay();
  } else {
    _dtheta.resize(100);
    _ddTheta = (thetaMax-_thetaMin)/100;
    _theta->eval(_thetaMin, thetaMax, _dtheta);
  }
}

GammaCompoundDensityObj::~GammaCompoundDensityObj() {
  // pass...
}

void
GammaCompoundDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  _k->mark();
  _theta->mark();
}

void
GammaCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();

  // Apply affine transform
  Tmin = Tmin-_shift;
  Tmax = Tmax-_shift;

  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double k = _kMin;
    for (int l=0; l<_dk.size(); l++, k+=_ddK) {
      double theta = _thetaMin;
      for (int m=0; m<_dtheta.size(); m++, theta+=_ddTheta) {
        out(i) += _ddK*_ddTheta * _dk[l] * _dtheta[m] *
            std::exp((k-1)*std::log(t) - t/theta -std::lgamma(k) -k*std::log(theta));
      }
    }
  }
}


void
GammaCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();
  // Apply affine transform
  Tmin = Tmin-_shift;
  Tmax = Tmax-_shift;
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double k = _kMin;
    for (int l=0; l<_dk.size(); l++, k+=_ddK) {
      double theta = _thetaMin;
      for (int m=0; m<_dtheta.size(); m++, theta+=_ddTheta) {
        out(i) += _ddK*_ddTheta * _dk[l] * _dtheta[m] *
            stochbb::gamma_li(k, t/theta);
      }
    }
  }
}

Density
GammaCompoundDensityObj::affine(double scale, double shift) const {
  return new GammaCompoundDensityObj(_k, *_theta->affine(scale, 0), scale*_shift+shift);
}

void
GammaCompoundDensityObj::rangeEst(double alpha, double &a, double &b) const {
  // Get quantiles for parameter distributions
  double a_k, b_k, a_theta, b_theta;
  _k->rangeEst(alpha, a_k, b_k);
  _theta->rangeEst(alpha, a_theta, b_theta);
  // Get maximum quantiles
  a = qgamma(alpha, a_k, a_theta);
  b = qgamma(1-alpha, b_k, b_theta);
  // Apply affine trafo on a & b;
  a = a+_shift;
  b = b+_shift;
}


/* ********************************************************************************************* *
 * Implementation of GammaCompoundObj
 * ********************************************************************************************* */
GammaCompoundObj::GammaCompoundObj(const Var &k, const Var &theta, const std::string &name)
  : CompoundObj(std::vector<Var> {k, theta}, name), _density(0)
{
  // Test for independence
  if (! k.mutuallyIndep(theta)) {
    AssumptionError err;
    err << "Cannot construct gamma-compound parameters are not independent.";
    throw err;
  }

  _density = new GammaCompoundDensityObj(k, theta);
  _density->unref();
}

void
GammaCompoundObj::mark() {
  if (isMarked()) { return; }
  CompoundObj::mark();
  if (_density) { _density->mark(); }
}

Density
GammaCompoundObj::density() {
  _density->mark(); return _density;
}



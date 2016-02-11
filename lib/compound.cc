#include "compound.hh"
#include "operators.hh"
#include "math.hh"
#include "logger.hh"

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
  : DensityObj(), _mu(0), _sigma(0), _muMin(0), _ddMu(0), _dmu(),
    _sigMin(0), _ddSig(0), _dsigma()
{
  _mu = *mu->density();
  _sigma = *sigma->density();
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
  if (_mu) _mu->mark();
  if (_sigma) _sigma->mark();
}

void
NormalCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
#pragma omp for
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
#pragma omp for
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
  : DensityObj(), _k(0), _theta(0), _shift(shift), _kMin(0), _ddK(0), _dk(),
    _thetaMin(0), _ddTheta(0), _dtheta()
{
  _k = *k->density();
  _theta = *theta->density();
  // Test for independence
  if (! k.mutuallyIndep(theta)) {
    AssumptionError err;
    err << "Cannot construct gamma-compound parameters are not independent.";
    throw err;
  }
  std::stringstream buffer;
  buffer << "Create Gamma compound with k=" << k << ", theta=" << theta << ".";
  logDebug() << buffer.str();
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

  // Check if theta is a delta distribution
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
  if (_k) _k->mark();
  if (_theta) _theta->mark();
}

void
GammaCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();

  // Apply affine transform
  Tmin = Tmin-_shift;
  Tmax = Tmax-_shift;

  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
#pragma omp for
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double k = _kMin;
    for (int l=0; l<_dk.size(); l++, k+=_ddK) {
      double theta = _thetaMin;
      for (int m=0; m<_dtheta.size(); m++, theta+=_ddTheta) {
        out(i) += _ddK*_ddTheta * _dk[l] * _dtheta[m] * stochbb::dgamma(t, k, theta);
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
#pragma omp for
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double k = _kMin;
    for (int l=0; l<_dk.size(); l++, k+=_ddK) {
      double theta = _thetaMin;
      for (int m=0; m<_dtheta.size(); m++, theta+=_ddTheta) {
        out(i) += _ddK*_ddTheta * _dk[l] * _dtheta[m] * stochbb::pgamma(t, k, theta);
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


/* ********************************************************************************************* *
 * Implementation of InvGammaCompoundDensityObj
 * ********************************************************************************************* */
InvGammaCompoundDensityObj::InvGammaCompoundDensityObj(DensityObj *alpha, DensityObj *beta, double shift)
  : DensityObj(), _alpha(alpha), _beta(beta), _shift(shift), _alphaMin(0), _ddAlpha(0), _dalpha(),
    _betaMin(0), _ddBeta(0), _dbeta()
{
  _init_int();
}

InvGammaCompoundDensityObj::InvGammaCompoundDensityObj(const Var &alpha, const Var &beta, double shift) throw (AssumptionError)
  : DensityObj(), _alpha(0), _beta(0), _shift(shift), _alphaMin(0), _ddAlpha(0), _dalpha(),
    _betaMin(0), _ddBeta(0), _dbeta()
{
  _alpha = *alpha->density();
  _beta = *beta->density();
  // Test for independence
  if (! alpha.mutuallyIndep(beta)) {
    AssumptionError err;
    err << "Cannot construct inverse gamma-compound: Parameters are not independent.";
    throw err;
  }
  _init_int();
}

void
InvGammaCompoundDensityObj::_init_int() {
  double alphaMax, betaMax;
  // Get intergration range for alpha
  _alpha->rangeEst(0.01, _alphaMin, alphaMax);
  // Get intergration range for beta
  _beta->rangeEst(0.01, _betaMin, betaMax);

  // Check if alpha is a delta distribution
  if (DeltaDensityObj *delta_alpha = dynamic_cast<DeltaDensityObj *>(_alpha)) {
    _dalpha.resize(1); _dalpha[0] = 1; _ddAlpha = 1; _alphaMin = delta_alpha->delay();
  } else {
    _dalpha.resize(100);
    _ddAlpha = (alphaMax-_alphaMin)/100;
    _alpha->eval(_alphaMin, alphaMax, _dalpha);
  }

  // Check if beta is a delta distribution
  if (DeltaDensityObj *delta_beta = dynamic_cast<DeltaDensityObj *>(_beta)) {
    _dbeta.resize(1); _dbeta[0] = 1; _ddBeta = 1; _betaMin = delta_beta->delay();
  } else {
    _dbeta.resize(100);
    _ddBeta = (betaMax-_betaMin)/100;
    _beta->eval(_betaMin, betaMax, _dbeta);
  }
}

InvGammaCompoundDensityObj::~InvGammaCompoundDensityObj() {
  // pass...
}

void
InvGammaCompoundDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  if (_alpha) _alpha->mark();
  if (_beta) _beta->mark();
}

void
InvGammaCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();

  // Apply affine transform
  Tmin = Tmin-_shift;
  Tmax = Tmax-_shift;

  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
#pragma omp for
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double alpha = _alphaMin;
    for (int l=0; l<_dalpha.size(); l++, alpha+=_ddAlpha) {
      double beta = _betaMin;
      for (int m=0; m<_dbeta.size(); m++, beta+=_ddBeta) {
        out(i) += _ddAlpha*_ddBeta * _dalpha[l] * _dbeta[m] * stochbb::dinvgamma(t, alpha, beta);
      }
    }
  }
}

void
InvGammaCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();
  // Apply affine transform
  Tmin = Tmin-_shift;
  Tmax = Tmax-_shift;
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
#pragma omp for
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double alpha = _alphaMin;
    for (int l=0; l<_dalpha.size(); l++, alpha+=_ddAlpha) {
      double beta = _betaMin;
      for (int m=0; m<_dbeta.size(); m++, beta+=_ddBeta) {
        out(i) += _ddAlpha*_ddBeta * _dalpha[l] * _dbeta[m] * stochbb::pinvgamma(t, alpha, beta);
      }
    }
  }
}

Density
InvGammaCompoundDensityObj::affine(double scale, double shift) const {
  return new InvGammaCompoundDensityObj(_alpha, *_beta->affine(scale, 0), scale*_shift+shift);
}

void
InvGammaCompoundDensityObj::rangeEst(double alpha, double &a, double &b) const {
  // Get quantiles for parameter distributions
  double a_alpha, b_alpha, a_beta, b_beta;
  _alpha->rangeEst(alpha, a_alpha, b_alpha);
  _beta->rangeEst(alpha, a_beta, b_beta);
  // Get maximum quantiles
  a = qinvgamma(alpha, a_alpha, a_beta);
  b = qinvgamma(1-alpha, b_alpha, b_beta);
  // Apply affine trafo on a & b;
  a = a+_shift;
  b = b+_shift;
}


/* ********************************************************************************************* *
 * Implementation of InvGammaCompoundObj
 * ********************************************************************************************* */
InvGammaCompoundObj::InvGammaCompoundObj(const Var &alpha, const Var &beta, const std::string &name)
  : CompoundObj(std::vector<Var> {alpha, beta}, name), _density(0)
{
  // Test for independence
  if (! alpha.mutuallyIndep(beta)) {
    AssumptionError err;
    err << "Cannot construct inverse gamma-compound parameters are not independent.";
    throw err;
  }

  _density = new InvGammaCompoundDensityObj(alpha, beta);
  _density->unref();
}

void
InvGammaCompoundObj::mark() {
  if (isMarked()) { return; }
  CompoundObj::mark();
  if (_density) { _density->mark(); }
}

Density
InvGammaCompoundObj::density() {
  _density->mark(); return _density;
}


/* ********************************************************************************************* *
 * Implementation of WeibullCompoundDensityObj
 * ********************************************************************************************* */
WeibullCompoundDensityObj::WeibullCompoundDensityObj(DensityObj *k, DensityObj *lambda, double shift)
  : DensityObj(), _k(k), _lambda(lambda), _shift(shift), _kMin(0), _ddK(0), _dk(),
    _lambdaMin(0), _ddLambda(0), _dlambda()
{
  _init_int();
}

WeibullCompoundDensityObj::WeibullCompoundDensityObj(const Var &k, const Var &lambda, double shift) throw (AssumptionError)
  : DensityObj(), _k(0), _lambda(0), _shift(shift), _kMin(0), _ddK(0), _dk(),
    _lambdaMin(0), _ddLambda(0), _dlambda()
{
  _k = *k->density();
  _lambda = *lambda->density();
  // Test for independence
  if (! k.mutuallyIndep(lambda)) {
    AssumptionError err;
    err << "Cannot construct Weibull-compound parameters are not independent.";
    throw err;
  }

  _init_int();
}

void
WeibullCompoundDensityObj::_init_int() {
  double kMax, lambdaMax;
  // Get intergration range for k
  _k->rangeEst(0.01, _kMin, kMax);
  // Get intergration range for theta
  _lambda->rangeEst(0.01, _lambdaMin, lambdaMax);

  // Check if theta is a delta distribution
  if (DeltaDensityObj *delta_k = dynamic_cast<DeltaDensityObj *>(_k)) {
    _dk.resize(1); _dk[0] = 1; _ddK = 1; _kMin = delta_k->delay();
  } else {
    _dk.resize(100);
    _ddK = (kMax-_kMin)/100;
    _k->eval(_kMin, kMax, _dk);
  }

  // Check if lambda is a delta distribution
  if (DeltaDensityObj *delta_lambda = dynamic_cast<DeltaDensityObj *>(_lambda)) {
    _dlambda.resize(1); _dlambda[0] = 1; _ddLambda = 1; _lambdaMin = delta_lambda->delay();
  } else {
    _dlambda.resize(100);
    _ddLambda = (lambdaMax-_lambdaMin)/100;
    _lambda->eval(_lambdaMin, lambdaMax, _dlambda);
  }
}

WeibullCompoundDensityObj::~WeibullCompoundDensityObj() {
  // pass...
}

void
WeibullCompoundDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  if (_k) _k->mark();
  if (_lambda) _lambda->mark();
}

void
WeibullCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();

  // Apply affine transform
  Tmin = Tmin-_shift;
  Tmax = Tmax-_shift;

  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
#pragma omp for
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double k = _kMin;
    for (int l=0; l<_dk.size(); l++, k+=_ddK) {
      double lambda = _lambdaMin;
      for (int m=0; m<_dlambda.size(); m++, lambda+=_ddLambda) {
        out(i) += _ddK*_ddLambda * _dk[l] * _dlambda[m] * stochbb::dweibull(t, k, lambda);
      }
    }
  }
}

void
WeibullCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  out.setZero();
  // Apply affine transform
  Tmin = Tmin-_shift;
  Tmax = Tmax-_shift;
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
#pragma omp for
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double k = _kMin;
    for (int l=0; l<_dk.size(); l++, k+=_ddK) {
      double lambda = _lambdaMin;
      for (int m=0; m<_dlambda.size(); m++, lambda+=_ddLambda) {
        out(i) += _ddK*_ddLambda * _dk[l] * _dlambda[m] * stochbb::pweibull(t, k, lambda);
      }
    }
  }
}

Density
WeibullCompoundDensityObj::affine(double scale, double shift) const {
  return new WeibullCompoundDensityObj(_k, *_lambda->affine(scale, 0), scale*_shift+shift);
}

void
WeibullCompoundDensityObj::rangeEst(double alpha, double &a, double &b) const {
  // Get quantiles for parameter distributions
  double a_k, b_k, a_lambda, b_lambda;
  _k->rangeEst(alpha, a_k, b_k);
  _lambda->rangeEst(alpha, a_lambda, b_lambda);
  // Get maximum quantiles
  a = qweibull(alpha, a_k, a_lambda);
  b = qweibull(1-alpha, b_k, b_lambda);
  // Apply affine trafo on a & b;
  a = a+_shift;
  b = b+_shift;
}


/* ********************************************************************************************* *
 * Implementation of WeibullCompoundObj
 * ********************************************************************************************* */
WeibullCompoundObj::WeibullCompoundObj(const Var &k, const Var &lambda, const std::string &name)
  : CompoundObj(std::vector<Var> {k, lambda}, name), _density(0)
{
  // Test for independence
  if (! k.mutuallyIndep(lambda)) {
    AssumptionError err;
    err << "Cannot construct Weibull-compound parameters are not independent.";
    throw err;
  }

  _density = new WeibullCompoundDensityObj(k, lambda);
  _density->unref();
}

void
WeibullCompoundObj::mark() {
  if (isMarked()) { return; }
  CompoundObj::mark();
  if (_density) { _density->mark(); }
}

Density
WeibullCompoundObj::density() {
  _density->ref(); return _density;
}



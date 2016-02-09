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
NormalCompoundDensityObj::NormalCompoundDensityObj(DensityObj *mu, DensityObj *sigma, double scale, double shift)
  : DensityObj(), _mu(mu), _sigma(sigma), _scale(scale), _shift(shift)
{
  // pass...
}

NormalCompoundDensityObj::NormalCompoundDensityObj(const Var &mu, const Var &sigma, double scale, double shift)
  : DensityObj(), _mu(*mu->density()), _sigma(*sigma->density()), _scale(scale), _shift(shift)
{
  // Test for independence
  if (! mu.mutuallyIndep(sigma)) {
    AssumptionError err;
    err << "Cannot construct normal-compound parameters are not independent.";
    throw err;
  }
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
  // Eval PDFs of mu & sigma
  const size_t steps = 100;
  Eigen::VectorXd dmu(steps), dsigma(steps);
  double muMin, muMax, sigMin, sigMax;
  _mu->rangeEst(0.01, muMin, muMax);
  _sigma->rangeEst(0.01, sigMin, sigMax);
  double ddmu = (muMax-muMin)/steps, ddsig = (sigMax-sigMin)/steps;
  _mu->eval(muMin, muMax, dmu);
  _sigma->eval(sigMin, sigMax, dsigma);

  out.setZero();
  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double mu = muMin;
    for (int k=0; k<steps; k++, mu+=ddmu) {
      double sigma = sigMin;
      for (int l=0; l<steps; l++, sigma+=ddsig) {
        out(i) += ddmu*ddsig * dmu[k] * dsigma[l] *
            std::exp( -(t-mu)*(t-mu)/(2*sigma*sigma)) / (std::sqrt(2*M_PI)*sigma*_scale);
      }
    }
  }
}


void
NormalCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Eval PDFs of mu & sigma
  const size_t steps = 100;
  Eigen::VectorXd dmu(steps), dsigma(steps);
  double muMin, muMax, sigMin, sigMax;
  _mu->rangeEst(0.01, muMin, muMax);
  _sigma->rangeEst(0.01, sigMin, sigMax);
  double ddmu = (muMax-muMin)/steps, ddsig = (sigMax-sigMin)/steps;
  _mu->eval(muMin, muMax, dmu);
  _sigma->eval(sigMin, sigMax, dsigma);

  out.setZero();
  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double mu = muMin;
    for (int k=0; k<steps; k++, mu+=ddmu) {
      double sigma = sigMin;
      for (int l=0; l<out.size(); l++, sigma+=ddsig) {
        out(i) += ddmu*dmu[k] * ddsig*dsigma[l] *
            0.5*(1+std::erf((t-mu)/(sigma*std::sqrt(2))));
      }
    }
  }
}

Density
NormalCompoundDensityObj::affine(double scale, double shift) const {
  return new NormalCompoundDensityObj(_mu, _sigma, scale*_scale, scale*_shift+shift);
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

  // Apply affine trafo on a & b;
  a = (a-_shift)/_scale;
  b = (b-_shift)/_scale;
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
GammaCompoundDensityObj::GammaCompoundDensityObj(DensityObj *k, DensityObj *theta, double scale, double shift)
  : DensityObj(), _k(k), _theta(theta), _scale(scale), _shift(shift)
{
  // pass...
}

GammaCompoundDensityObj::GammaCompoundDensityObj(const Var &k, const Var &theta, double scale, double shift)
  : DensityObj(), _k(*k->density()), _theta(*theta->density()), _scale(scale), _shift(shift)
{
  // Test for independence
  if (! k.mutuallyIndep(theta)) {
    AssumptionError err;
    err << "Cannot construct gamma-compound parameters are not independent.";
    throw err;
  }
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
  // Eval PDFs of k & theta
  Eigen::VectorXd dk(out.size()), dtheta(out.size());
  _k->eval(Tmin, Tmax, dk); _theta->eval(Tmin, Tmax, dtheta);

  out.setZero();
  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double k = Tmin;
    for (int l=0; l<out.size(); l++, k+=dt) {
      double theta = Tmin;
      for (int m=0; m<out.size(); m++, theta+=dt) {
        out(i) += dt*dt * dk[l] * dtheta[m] *
            std::exp((k-1)*std::log(t) - t/theta -std::lgamma(k) -k*std::log(theta*_scale));
      }
    }
  }
}


void
GammaCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Eval PDFs of k & theta
  Eigen::VectorXd dk(out.size()), dtheta(out.size());
  _k->eval(Tmin, Tmax, dk); _theta->eval(Tmin, Tmax, dtheta);

  out.setZero();
  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double k = Tmin;
    for (int l=0; l<out.size(); l++, k+=dt) {
      double theta = Tmin;
      for (int m=0; m<out.size(); m++, theta+=dt) {
        out(i) += dt*dt * dk[l] * dtheta[m] *
            stochbb::gamma_li(k, t/theta);
      }
    }
  }
}

Density
GammaCompoundDensityObj::affine(double scale, double shift) const {
  return new GammaCompoundDensityObj(_k, _theta, scale*_scale, scale*_shift+shift);
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
  a = (a-_shift)/_scale;
  b = (b-_shift)/_scale;
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



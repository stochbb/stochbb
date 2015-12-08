#include "compound.hh"
#include "operators.hh"
#include "math.hh"


using namespace sbb;

// Helper function to construct a vector of variables inline
inline std::vector<Var> mk_var_vec(const Var &a, const Var &b) {
  std::vector<Var> vec; vec.reserve(2);
  vec.push_back(a); vec.push_back(b);
  return vec;
}


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
NormalCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  // Eval PDFs of mu & sigma
  Eigen::VectorXd dmu(out.size()), dsigma(out.size());
  _mu->eval(Tmin, Tmax, dmu); _sigma->eval(Tmin, Tmax, dsigma);

  out.setZero();
  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double mu = Tmin;
    for (int k=0; k<out.size(); k++, mu+=dt) {
      double sigma = Tmin;
      for (int l=0; l<out.size(); l++, sigma+=dt) {
        out(i) += dt*dt * dmu[k] * dsigma[l] *
            std::exp( -(t-mu)*(t-mu)/(2*sigma*sigma)) / (std::sqrt(2*M_PI)*sigma);
      }
    }
  }
}


void
NormalCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  // Eval PDFs of mu & sigma
  Eigen::VectorXd dmu(out.size()), dsigma(out.size());
  _mu->eval(Tmin, Tmax, dmu); _sigma->eval(Tmin, Tmax, dsigma);

  out.setZero();
  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;
  // For each time t
  double t=Tmin, dt = (Tmax-Tmin)/out.size();
  for (int i=0; i<out.size(); i++, t+=dt) {
    // Perform numerical integral over coeffs
    double mu = Tmin;
    for (int k=0; k<out.size(); k++, mu+=dt) {
      double sigma = Tmin;
      for (int l=0; l<out.size(); l++, sigma+=dt) {
        out(i) += dt*dmu[k] * dt*dsigma[l] *
            0.5*(1+std::erf((t-mu)/(sigma*std::sqrt(2))));
      }
    }
  }
}

Density
NormalCompoundDensityObj::affine(double scale, double shift) const {
  return new NormalCompoundDensityObj(_mu, _sigma, scale*_scale, scale*_shift+shift);
}


/* ********************************************************************************************* *
 * Implementation of NormalCompoundObj
 * ********************************************************************************************* */
NormalCompoundObj::NormalCompoundObj(const Var &mu, const Var &sigma, const std::string &name)
  : CompoundObj(mk_var_vec(mu, sigma), name), _density(0)
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
GammaCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
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
            std::exp((k-1)*std::log(t) - t/theta -std::lgamma(k) -k*std::log(theta));
      }
    }
  }
}


void
GammaCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
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
            sbb::gamma_li(k, t/theta);
      }
    }
  }
}

Density
GammaCompoundDensityObj::affine(double scale, double shift) const {
  return new GammaCompoundDensityObj(_k, _theta, scale*_scale, scale*_shift+shift);
}


/* ********************************************************************************************* *
 * Implementation of GammaCompoundObj
 * ********************************************************************************************* */
GammaCompoundObj::GammaCompoundObj(const Var &k, const Var &theta, const std::string &name)
  : CompoundObj(mk_var_vec(k, theta), name), _density(0)
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



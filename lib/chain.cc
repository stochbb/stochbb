#include "chain.hh"
#include "exception.hh"
#include <unsupported/Eigen/FFT>
#include <list>
#include <algorithm>
#include <complex>
#include "operators.hh"
#include "logger.hh"

using namespace stochbb;

// tiny helper function to sort vectors of densities
inline int density_compare(const Density &a, const Density &b) {
  return a.compare(b);
}

/* ********************************************************************************************* *
 * combine densities
 * ********************************************************************************************* */
// Checks if two densities can be combined
bool convolution_can_combine(const Density &a, const Density &b) {
  // Check by type
  if (dynamic_cast<DeltaDensityObj *>(*a)) {
    // If LHS is a delta distribution -> yes
    return true;
  } else if (dynamic_cast<UniformDensityObj *>(*a)) {
    // If LHS is a uniform density ...
    if (dynamic_cast<DeltaDensityObj *>(*b)) {
      // ... and RHS is a delta density -> yes
      return true;
    }
  } else if (dynamic_cast<NormalDensityObj *>(*a)) {
    // If RHS is a normal distribution ...
    if (dynamic_cast<DeltaDensityObj *>(*b)) {
      // ... and RHS is a delta density -> yes
      return true;
    } else if (dynamic_cast<NormalDensityObj *>(*b)) {
      // ... and RHS is a normal density -> yes
      return true;
    }
  } else if (GammaDensityObj *gamma_a = dynamic_cast<GammaDensityObj *>(*a)) {
    // If LHS is a gamma density ...
    if (dynamic_cast<DeltaDensityObj *>(*b)) {
      // If RHS is a delta density ...
      return true;
    } else if (GammaDensityObj *gamma_b = dynamic_cast<GammaDensityObj *>(*b)) {
      // If LHS is a gamma too and has the same scale
      return gamma_a->theta() == gamma_b->theta();
    }
  }
  // Otherwise -> no
  return false;
}

// Derives new the analytic convolution of densities (if possible, check with
// convolution_can_combine). Returns a new reference to a convolution object
Density
convolution_combine(const Density &a, const Density &b) {
  logDebug() << "Convolve densities...";
  if (DeltaDensityObj *delta_a = dynamic_cast<DeltaDensityObj *>(*a)) {
    // If LHS is delta -> turn into affine trafo
    return b.affine(1, delta_a->delay());
  } else if (UniformDensityObj *unif_a = dynamic_cast<UniformDensityObj *>(*a)) {
    // If LHS is a uniform density ...
    if (DeltaDensityObj *delta_b = dynamic_cast<DeltaDensityObj *>(*b)) {
      // ... and RHS is a delta density -> yes
      return new UniformDensityObj(unif_a->a()+delta_b->delay(),
                                   unif_a->b()+delta_b->delay());
    }
  } else if (NormalDensityObj *norm_a = dynamic_cast<NormalDensityObj *>(*a)) {
    // If LHS is a normal distribution ...
    if (DeltaDensityObj *delta_b = dynamic_cast<DeltaDensityObj *>(*b)) {
      // ... and RHS is a delta density
      return new NormalDensityObj(norm_a->mu()+delta_b->delay(), norm_a->sigma());
    } else if (NormalDensityObj *norm_b = dynamic_cast<NormalDensityObj *>(*b)) {
      // ... and RHS is a normal density too
      return new NormalDensityObj(norm_a->mu()+norm_b->mu(),
                                  std::sqrt(norm_a->sigma()*norm_a->sigma()
                                            + norm_b->sigma()*norm_b->sigma()));
    }
  } else if (GammaDensityObj *gamma_a = dynamic_cast<GammaDensityObj *>(*a)) {
    // If LHS is a gamma density ...
    if (DeltaDensityObj *delta_b = dynamic_cast<DeltaDensityObj *>(*b)) {
      // ... and RHS is a delta density
      return new GammaDensityObj(gamma_a->k(), gamma_a->theta(), gamma_a->shift()+delta_b->delay());
    } else if (GammaDensityObj *gamma_b = dynamic_cast<GammaDensityObj *>(*b)) {
      // If LHS is a gamma too and has the same scale
      if (gamma_a->theta() == gamma_b->theta()) {
        return new GammaDensityObj(gamma_a->k()+gamma_b->k(), gamma_a->theta(),
                                   gamma_a->shift()+gamma_b->shift());
      }
    }
  }
  // Otherwise --> oops
  AssumptionError err;
  err << "Cannot combine densities.";
  throw err;
}

Density
stochbb::convolve(const std::vector<Density> &densities, double scale, double shift) {
  // copy vector
  std::vector<Density> dens(densities);
  // Sort densities w.r.t type and parameters
  std::sort(dens.begin(), dens.end(), density_compare);

  logDebug() << "Try to convolve " << dens.size() << " densities...";
  // Try to combine some of the densities
  std::vector<Density>::iterator last = dens.begin();
  std::vector<Density>::iterator current = dens.begin(); current++;
  while (current != dens.end()) {
    if (convolution_can_combine(*last, *current)) {
      // If densities can be combined -> combine & replace last density
      *last = convolution_combine(*last, *current);
      // erase combined density
      current = dens.erase(current);
    } else {
      // If densities cannot be combined -> advance iterators
      last++; current++;
    }
  }
  logDebug() << "... done.";

  if (1 == dens.size()) {
    // If only one density is left -> unpack
    return dens.back();
  }

  // Otherwise construct convolution density from list of densities
  return new ConvolutionDensityObj(dens);
}


/* ********************************************************************************************* *
 * Implementation of ConvolutionDensityObj
 * ********************************************************************************************* */
ConvolutionDensityObj::ConvolutionDensityObj(const std::vector<DensityObj *> &densities, double scale, double shift)
  : DensityObj(), _densities(densities), _scale(scale), _shift(shift)
{
  // pass...
}

ConvolutionDensityObj::ConvolutionDensityObj(const std::vector<Density> &densities, double scale, double shift)
  : DensityObj(), _densities(), _scale(scale), _shift(shift)
{
  _densities.reserve(densities.size());
  for (size_t i=0; i<densities.size(); i++) {
    _densities.push_back(*densities[i]);
  }
}

ConvolutionDensityObj::ConvolutionDensityObj(const std::vector<Var> &variables, double scale, double shift)
  : DensityObj(), _densities(), _scale(scale), _shift(shift)
{
  // Get & store the densities of all variables, assuming they are mutually independent.
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(*variables[i]->density());
  }
}

ConvolutionDensityObj::~ConvolutionDensityObj() {
  // pass...
}

void
ConvolutionDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->mark();
  }
}

void
ConvolutionDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;
  // Get sample period
  double dt = (Tmax-Tmin)/out.size();
  // Get sample rate
  double Fs = out.size()/(Tmax-Tmin);
  // Compute time-shift
  std::complex<double> dphi(0, 2*M_PI*Tmin*Fs/double(2*out.size()));

  // Allocate some buffers & FFT trafo
  Eigen::FFT<double> fft;
  Eigen::VectorXd tmp1(2*out.size());  tmp1.setZero();
  Eigen::VectorXcd tmp2(2*out.size());
  Eigen::VectorXcd prod(2*out.size()); prod.setOnes();

  // Perform FFT convolution:
  // For each density ...
  for (size_t i=0; i<_densities.size(); i++) {
    // ... eval
    _densities[i]->eval(Tmin, Tmax, out);
    // ... turn into PMF
    tmp1.head(out.size()) = out*dt;
    // ... tmp2 = FFT(tmp1)
    fft.fwd(tmp2, tmp1);
    // Skip first PDF, this avoids the back-shift of the product
    if (i>0) {
      // apply time shift
      for (int j=1; j<out.size();j++) {
        tmp2[j] *= std::exp(-double(j)*dphi);
        tmp2[2*out.size()-j] *= std::exp(double(j)*dphi);
      }
    }
    // prod = prod * FFT( density )
    prod.array() *= tmp2.array();
  }
  // tmp1 = InvFFT(prod)
  fft.inv(tmp1, prod);
  // Store result
  out.noalias() = tmp1.head(out.size())/(dt*_scale);
}

void
ConvolutionDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  if (0 == out.size()) { return; }
  // Get sample period
  double dt = (Tmax-Tmin)/out.size();
  // Eval PDF
  this->eval(Tmin, Tmax, out);
  // Rescale first element
  out[0] *= dt;
  // compute cummulative
  for (int i=1; i<out.size(); i++) {
    out[i] = out[i-1] + out[i]*dt;
  }
}

Density
ConvolutionDensityObj::affine(double scale, double shift) const {
  return new ConvolutionDensityObj(_densities, _scale*scale, scale*_shift+shift);
}

void
ConvolutionDensityObj::rangeEst(double alpha, double &a, double &b) const {
  /* This function ties to approximate the alpha-quantiles of the convolution
   * of the densities. The approx. quantiles are obtained as the support of the convolution
   * of characteristic functions of the alpha-quantiles of the distributions. This can be
   * a relatively bad approximation. */

  // Get quantiles of the first density
  _densities[0]->rangeEst(alpha,a,b);
  double amin=a, bmax=b;
  for (size_t i=1; i<_densities.size(); i++) {
    // Get quantiles of next desity
    double c,d; _densities[i]->rangeEst(alpha, c,d);
    // Update quantiles;
    a += c; b += d;
    amin = std::min(amin, c);
    bmax = std::max(bmax, d);
  }
  // apply scale & shift on a & b:
  a = _scale*a+_shift;
  b = _scale*b+_shift;

  // Finally update amin bmax with a,b
  a = std::min(a, amin);
  b = std::max(b, bmax);
}

int
ConvolutionDensityObj::compare(const DensityObj &other) const {
  // Compare by type
  if (int res = DensityObj::compare(other)) { return res; }
  // If types match
  const ConvolutionDensityObj *o_conv = dynamic_cast<const ConvolutionDensityObj *>(&other);
  // Compare by number of densities
  if (_densities.size() < o_conv->_densities.size()) { return -1; }
  else if (_densities.size() > o_conv->_densities.size()) { return 1; }
  // Compare densities
  for (size_t i=0; i<_densities.size(); i++) {
    if (int res = _densities[i]->compare(*o_conv->_densities[i])) { return res; }
  }
  // equal
  return 0;
}

void
ConvolutionDensityObj::print(std::ostream &stream) const {
  stream << "<ConvolutionDensityObj of";
  for (size_t i=0; i<_densities.size(); i++) {
    stream << " "; _densities[i]->print(stream);
  }
  stream << " shift=" << _shift
         << " scale=" << _scale
         << " #" << this << ">";
}


/* ********************************************************************************************* *
 * Implementation of ChainObj
 * ********************************************************************************************* */
ChainObj::ChainObj(const std::vector<Var> &variables, const std::string &name)
  : DerivedVarObj(variables, name), _density(0)
{
  // Check for independence
  if (! independent(variables)) {
    AssumptionError err;
    err << "Cannot create chain: Variables not independent.";
    throw err;
  }

  // Get vector of densities
  std::vector<Density> dens; dens.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    dens.push_back(variables[i].density());
  }

  // Assemble convolution density
  _density = *convolve(dens);
  //_density->unref();
}

ChainObj::~ChainObj() {
  // pass...
}

void
ChainObj::mark() {
  if (isMarked()) { return; }
  DerivedVarObj::mark();
  // mark density
  if (_density) {
    _density->mark();
  }
}

Density ChainObj::density() {
  _density->ref();
  return _density;
}

void
ChainObj::print(std::ostream &stream) const {
  stream << "<ChainObj";
  for (size_t i=0; i<_variables.size(); i++) {
    stream << " "; _variables[i]->print(stream);
  }
  stream << " #" << this << ">";
}

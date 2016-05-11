#include "mixture.hh"
#include "exception.hh"
#include "rng.hh"

using namespace stochbb;

inline size_t
_find_index(double p, size_t a, size_t b, const std::vector<double> &cdf) {
  while (1 < (b-a)) {
    size_t mid = a+(b-a)/2;
    if (p < cdf[mid]) { b = mid; }
    else { a = mid; }
  }
  return b;
}


/* ********************************************************************************************* *
 * Implementation of MixtrueDensityObj
 * ********************************************************************************************* */
MixtureDensityObj::MixtureDensityObj(const std::vector<double> &weights, const std::vector<DensityObj *> &densities,
                                     double scale, double shift)
  : DensityObj(), _weights(weights), _densities(densities), _scale(scale), _shift(shift)
{
  // pass...
}

MixtureDensityObj::MixtureDensityObj(const std::vector<double> &weights, const std::vector<VarObj *> &variables,
                                     double scale, double shift)
  : DensityObj(), _weights(weights), _densities(), _scale(scale), _shift(shift)
{
  // Normalize weights
  double sum = 0;
  for (size_t i=0; i<_weights.size(); i++) { sum += _weights[i]; }
  for (size_t i=0; i<_weights.size(); i++) { _weights[i] /= sum; }
  // Store densities
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(*variables[i]->density());
  }
}

void
MixtureDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->mark();
  }
}

void
MixtureDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // apply affine transform
  Tmin = (Tmin - _shift)/_scale;
  Tmax = (Tmax - _shift)/_scale;

  out.setZero();
  Eigen::VectorXd tmp(out.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->eval(Tmin, Tmax, tmp);
    out += _weights[i]*tmp/_scale;
  }
}

void
MixtureDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // apply affine transform
  Tmin = (Tmin - _shift)/_scale;
  Tmax = (Tmax - _shift)/_scale;

  out.setZero();
  Eigen::VectorXd tmp(out.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->evalCDF(Tmin, Tmax, tmp);
    out += _weights[i]*tmp;
  }
}

Density
MixtureDensityObj::affine(double scale, double shift) const {
  return new MixtureDensityObj(_weights, _densities, scale*_scale, scale*_shift + shift);
}

void
MixtureDensityObj::rangeEst(double alpha, double &a, double &b) const {
  // Returns an estimate of the quantile as the union of all quantiles
  _densities[0]->rangeEst(alpha, a, b);
  for (size_t i=0; i<_densities.size(); i++) {
    double c,d; _densities[i]->rangeEst(alpha, c, d);
    a = std::min(a, c);
    b = std::max(b, d);
  }
}

int
MixtureDensityObj::compare(const DensityObj &other) const {
  // compare by type
  if (int cmp = DensityObj::compare(other))
    return cmp;
  const MixtureDensityObj &oobj = dynamic_cast<const MixtureDensityObj &>(other);

  // compare by number of mixed densities.
  if (_densities.size() < oobj._densities.size()) { return -1; }
  else if (_densities.size() > oobj._densities.size()) { return 1; }

  // compare mixed densities element by element
  for (size_t i=0; i<_densities.size(); i++) {
    if (int cmp = _densities[i]->compare(* oobj._densities[i]))
      return cmp;
    if (_weights[i] < oobj._weights[i])
      return -1;
    if (_weights[i] > oobj._weights[i])
      return 1;
  }

  return 0;
}

void
MixtureDensityObj::print(std::ostream &stream) const {
  stream << "<MixtureDensity of";
  for (size_t i=0; i<_densities.size(); i++) {
    stream << " " << _weights[i] << "*"; _densities[i]->print(stream);
  }
  if (_shift)
    stream << " shift=" << _shift;
  if (1 != _scale)
    stream << " scale=" << _scale;
  stream << " #" << this << ">";
}


/* ********************************************************************************************* *
 * Implementation of MixtrueObj
 * ********************************************************************************************* */
MixtureObj::MixtureObj(const std::vector<double> &weights, const std::vector<Var> &variables,
                       const std::string &name) throw (Error)
  : DerivedVarObj(variables, name), _weights(weights), _density(0)
{
  if (_weights.size() != variables.size()) {
    AssumptionError err;
    err << "Cannot create MixtureObj, number of weights must match number of variables.";
    throw err;
  }

  _density = new MixtureDensityObj(_weights, _variables);
  _density->unref();
}

void
MixtureObj::mark() {
  if (isMarked()) { return; }
  DerivedVarObj::mark();
  if (_density) { _density->mark(); }
}

Density
MixtureObj::density() {
  _density->ref();
  return _density;
}


void
MixtureObj::sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                   Eigen::Ref<Eigen::MatrixXd> samples) const
{
  // compute cummulative distribution of weights
  std::vector<double> cum; cum.reserve(this->numVariables());
  for (size_t i=0; i<this->numVariables(); i++) {
    cum[i] = (i>0) ? (cum[i-1]+this->weight(i)) : this->weight(i);
  }
  // normalize cdf
  for (size_t i=0; i<this->numVariables(); i++) {
    cum[i] /= cum[this->numVariables()-1];
  }
  RNG &rng = RNG::get();
  std::uniform_real_distribution<double> sampler(0,1);
  for (int i=0; i<samples.rows(); i++) {
    // select a variable randomly
    double p = sampler(rng);
    size_t idx = (p < cum[0]) ? 0 : _find_index(p, 0, this->numVariables()-1, cum);
    // select sample
    samples(i, outIdx) = samples(i, indices(idx));
  }
}

void
MixtureObj::print(std::ostream &stream) const {
  stream << "<Mixture of";
  for (size_t i=0; i<_variables.size(); i++) {
    stream << " " << _weights[i] << "*"; _variables[i]->print(stream);
  }
  stream << " density="; _density->print(stream);
  stream << " #" << this << ">";
}


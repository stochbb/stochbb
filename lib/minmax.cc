#include "minmax.hh"
#include "exception.hh"
#include "chain.hh"

using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of MaximumDensityObj
 * ********************************************************************************************* */
MaximumDensityObj::MaximumDensityObj(const std::vector<DensityObj *> &densities, double scale, double shift)
  : DensityObj(), _densities(densities), _scale(scale), _shift(shift)
{
  // pass...
}

MaximumDensityObj::MaximumDensityObj(const std::vector<VarObj *> &variables, double scale, double shift)
  : DensityObj(), _densities(), _scale(scale), _shift(shift)
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(*variables[i]->density());
  }
}

MaximumDensityObj::MaximumDensityObj(const std::vector<Var> &variables, double scale, double shift)
  : DensityObj(), _densities(), _scale(scale), _shift(shift)
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(*variables[i].density());
  }
}

MaximumDensityObj::~MaximumDensityObj() {
  // pass...
}

void
MaximumDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->mark();
  }
}

void
MaximumDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;

  Eigen::VectorXd tmp(out.size());
  Eigen::MatrixXd pdfs(out.size(), _densities.size());
  Eigen::MatrixXd cdfs(out.size(), _densities.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->eval(Tmin, Tmax, tmp);    pdfs.col(i) = tmp;
    _densities[i]->evalCDF(Tmin, Tmax, tmp); cdfs.col(i) = tmp;
  }

  out.setZero();
  for (size_t i=0; i<_densities.size(); i++) {
    tmp.setOnes();
    for (size_t j=0; j<_densities.size(); j++) {
      if (i == j) {
        tmp.array() *= pdfs.col(j).array();
      } else {
        tmp.array() *= cdfs.col(j).array();
      }
    }
    out += tmp/_scale;
  }
}

void
MaximumDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  Eigen::VectorXd tmp(out.size());

  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;

  out.setOnes();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->evalCDF(Tmin, Tmax, tmp);
    out.array() *= tmp.array();
  }
}

Density
MaximumDensityObj::affine(double scale, double shift) const {
  return new MaximumDensityObj(_densities, _scale*scale, scale*_shift+shift);
}

void
MaximumDensityObj::rangeEst(double alpha, double &a, double &b) const {
  // Obtain an estimate of the alpha-quantiles as the maxiumum of all quantile intervals
  _densities[0]->rangeEst(alpha, a, b);
  for (size_t i=1; i<_densities.size(); i++) {
    double c,d; _densities[i]->rangeEst(alpha, c,d);
    a = std::max(a, c);
    b = std::max(b, d);
  }
}

int
MaximumDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // If types match
  const MaximumDensityObj *o_max = dynamic_cast<const MaximumDensityObj *>(&other);
  // Compare by number of densities
  if (_densities.size() < o_max->_densities.size()) { return -1; }
  else if (_densities.size() > o_max->densities().size()) { return 1; }
  // Compare densities
  for (size_t i=0; i<_densities.size(); i++) {
    if (int res = _densities[i]->compare(*o_max->_densities[i])) { return res; }
  }
  // otherwise -> identical
  return 0;
}

void
MaximumDensityObj::print(std::ostream &stream) const {
  stream << "<MaximumDensity of";
  for (size_t i=0; i<_densities.size(); i++) {
    stream << " "; _densities[i]->print(stream);
  }
  if (_shift)
    stream << " shift=" << _shift;
  if (1 != _scale)
    stream << " scale=" << _scale;
  stream << " #" << this << ">";
}


/* ********************************************************************************************* *
 * Implementation of MinimumDensityObj
 * ********************************************************************************************* */
MinimumDensityObj::MinimumDensityObj(const std::vector<DensityObj *> &densities, double scale, double shift)
  : DensityObj(), _densities(densities), _scale(scale), _shift(shift)
{
  // pass...
}

MinimumDensityObj::MinimumDensityObj(const std::vector<VarObj *> &variables, double scale, double shift)
  : DensityObj(), _densities(), _scale(scale), _shift(shift)
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(*variables[i]->density());
  }
}

MinimumDensityObj::MinimumDensityObj(const std::vector<Var> &variables, double scale, double shift)
  : DensityObj(), _densities(), _scale(scale), _shift(shift)
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(*variables[i].density());
  }
}

MinimumDensityObj::~MinimumDensityObj() {
  // pass...
}

void
MinimumDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->mark();
  }
}

void
MinimumDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;

  Eigen::MatrixXd pdfs(out.size(), _densities.size());
  Eigen::MatrixXd cdfs(out.size(), _densities.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->eval(Tmin, Tmax, pdfs.col(i));
    _densities[i]->evalCDF(Tmin, Tmax, cdfs.col(i));
  }

  Eigen::VectorXd tmp(out.size());
  out.setZero();
  for (size_t i=0; i<_densities.size(); i++) {
    tmp.setOnes();
    for (size_t j=0; j<_densities.size(); j++) {
      if (i == j) {
        tmp.array() *= (-pdfs.col(j).array());
      } else {
        tmp.array() *= (1-cdfs.col(j).array());
      }
    }
    out += tmp;
  }
  out = -out/_scale;
}

void
MinimumDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  Eigen::VectorXd tmp(out.size());

  // Apply affine transform
  Tmin = (Tmin-_shift)/_scale;
  Tmax = (Tmax-_shift)/_scale;

  out.setOnes();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->evalCDF(Tmin, Tmax, tmp);
    out.array() *= (1-tmp.array());
  }
  out.array() = (1-out.array());
}

Density
MinimumDensityObj::affine(double scale, double shift) const {
  return new MinimumDensityObj(_densities, scale*_scale, scale*_shift+shift);
}

void
MinimumDensityObj::rangeEst(double alpha, double &a, double &b) const {
  // Obtain an estimate of the alpha-quantiles as the minimum of all quantile intervals
  _densities[0]->rangeEst(alpha, a, b);
  for (size_t i=1; i<_densities.size(); i++) {
    double c,d; _densities[i]->rangeEst(alpha, c,d);
    a = std::min(a, c);
    b = std::min(b, d);
  }
}

int
MinimumDensityObj::compare(const DensityObj &other) const {
  // Compare types
  if (int res = DensityObj::compare(other)) { return res; }
  // If types match
  const MinimumDensityObj *o_max = dynamic_cast<const MinimumDensityObj *>(&other);
  // Compare by number of densities
  if (_densities.size() < o_max->_densities.size()) { return -1; }
  else if (_densities.size() > o_max->densities().size()) { return 1; }
  // Compare densities
  for (size_t i=0; i<_densities.size(); i++) {
    if (int res = _densities[i]->compare(*o_max->_densities[i])) { return res; }
  }
  // otherwise -> identical
  return 0;
}

void
MinimumDensityObj::print(std::ostream &stream) const {
  stream << "<MinimumDensity of";
  for (size_t i=0; i<_densities.size(); i++) {
    stream << " "; _densities[i]->print(stream);
  }
  if (_shift)
    stream << " shift=" << _shift;
  if (1 != _scale)
    stream << " scale=" << _scale;
  stream << " #" << this << ">";
}


/* ********************************************************************************************* *
 * Implementation of MaximumObj
 * ********************************************************************************************* */
MaximumObj::MaximumObj(const std::vector<Var> &variables, const std::string &name)
  : DerivedVarObj(variables, name), _density(0)
{
  // Construct density
  _density = new MaximumDensityObj(_variables);
  _density->unref();
}

MaximumObj::~MaximumObj() {
  // pass...
}

void
MaximumObj::mark() {
  if (isMarked()) { return; }
  VarObj::mark();
  if (_density) { _density->mark(); }
}

Density
MaximumObj::density() {
  _density->ref();
  return _density;
}

void
MaximumObj::sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                   Eigen::Ref<Eigen::MatrixXd> samples) const
{
  for (int i=0; i<samples.rows(); i++) {
    samples(i,outIdx) = samples(i, indices(0));
    for (size_t j=1; j<this->numVariables(); j++) {
      samples(i,outIdx) = std::max(samples(i,outIdx), samples(i, indices(j)));
    }
  }
}

void
MaximumObj::print(std::ostream &stream) const {
  stream << "<Maximum of";
  for (size_t i=0; i<_variables.size(); i++) {
    stream << " "; _variables[i]->print(stream);
  }
  stream << " density="; _density->print(stream);
  stream << " #" << this << ">";
}


/* ********************************************************************************************* *
 * Implementation of MinimumObj
 * ********************************************************************************************* */
MinimumObj::MinimumObj(const std::vector<Var> &variables, const std::string &name)
  : DerivedVarObj(variables, name), _density(0)
{
  // Construct density
  _density = new MinimumDensityObj(_variables);
  _density->unref();
}

MinimumObj::~MinimumObj() {
  // pass...
}

void
MinimumObj::mark() {
  if (isMarked()) { return; }
  VarObj::mark();
  if (_density) { _density->mark(); }
}

Density
MinimumObj::density() {
  if (_density)
    _density->ref();
  return _density;
}

void
MinimumObj::sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                   Eigen::Ref<Eigen::MatrixXd> samples) const
{
  for (int i=0; i<samples.rows(); i++) {
    samples(i,outIdx) = samples(i, indices(0));
    for (size_t j=1; j<this->numVariables(); j++) {
      samples(i,outIdx) = std::min(samples(i,outIdx), samples(i, indices(j)));
    }
  }
}

void
MinimumObj::print(std::ostream &stream) const {
  stream << "<Minimum of";
  for (size_t i=0; i<_variables.size(); i++) {
    stream << " "; _variables[i]->print(stream);
  }
  stream << " density="; _density->print(stream);
  stream << " #" << this << ">";
}

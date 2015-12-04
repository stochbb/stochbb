#include "minmax.hh"
#include "exception.hh"
#include "chain.hh"

using namespace sbb;


/* ********************************************************************************************* *
 * Implementation of MaximumDensityObj
 * ********************************************************************************************* */
MaximumDensityObj::MaximumDensityObj(const std::vector<VarObj *> &variables)
  : DensityObj(), _densities()
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(*variables[i]->density());
  }
}

MaximumDensityObj::MaximumDensityObj(const std::vector<Var> &variables)
  : DensityObj(), _densities()
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
MaximumDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
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
    out += tmp;
  }
}

void
MaximumDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  Eigen::VectorXd tmp(out.size());

  out.setOnes();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->evalCDF(Tmin, Tmax, tmp);
    out.array() *= tmp.array();
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
  stream << "<MaximumDensityObj of";
  for (size_t i=0; i<_densities.size(); i++) {
    stream << " "; _densities[i]->print(stream);
  }
  stream << " #" << this << ">";
}


/* ********************************************************************************************* *
 * Implementation of MinimumDensityObj
 * ********************************************************************************************* */
MinimumDensityObj::MinimumDensityObj(const std::vector<VarObj *> &variables)
  : DensityObj(), _densities()
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(*variables[i]->density());
  }
}

MinimumDensityObj::MinimumDensityObj(const std::vector<Var> &variables)
  : DensityObj(), _densities()
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
MinimumDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  double dx = (Tmax-Tmin)/out.size();
  Eigen::VectorXd tmp(out.size());
  Eigen::MatrixXd pdfs(out.size(), _densities.size());
  Eigen::MatrixXd cdfs(out.size(), _densities.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->eval(Tmin, Tmax, tmp);    pdfs.col(i) = tmp*dx;
    _densities[i]->evalCDF(Tmin, Tmax, tmp); cdfs.col(i) = tmp;
  }

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
  out = -out/dx;
}

void
MinimumDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  Eigen::VectorXd tmp(out.size());

  out.setOnes();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->evalCDF(Tmin, Tmax, tmp);
    out.array() *= (1-tmp.array());
  }
  out.array() = (1-out.array());
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
  stream << "<MinimumDensityObj of";
  for (size_t i=0; i<_densities.size(); i++) {
    stream << " "; _densities[i]->print(stream);
  }
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


/* ********************************************************************************************* *
 * Implementation of MinimumObj
 * ********************************************************************************************* */
MinimumObj::MinimumObj(const std::vector<Var> &variables, const std::string &name)
  : DerivedVarObj(variables, name), _density(0)
{
  // Construct density
  _density = new MinimumDensityObj(_variables);
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

Density MinimumObj::density() {
  _density->ref();
  return _density;
}

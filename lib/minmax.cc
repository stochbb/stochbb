#include "minmax.hh"
#include "exception.hh"

using namespace sbb;


/* ********************************************************************************************* *
 * Implementation of MaximumDensityObj
 * ********************************************************************************************* */
MaximumDensityObj::MaximumDensityObj(const std::vector<VarObj *> &variables)
  : DensityObj(), _densities()
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(variables[i]->density());
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

void
MaximumDensityObj::sample(Eigen::VectorXd &out) const {
  Eigen::VectorXd tmp;
  Eigen::MatrixXd samples(out.size(), _densities.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->sample(tmp); samples.col(i) = tmp;
  }
  for (size_t i=0; i<out.size(); i++) {
    out[i] = samples.row(i).maxCoeff();
  }
}


/* ********************************************************************************************* *
 * Implementation of MinimumDensityObj
 * ********************************************************************************************* */
MinimumDensityObj::MinimumDensityObj(const std::vector<VarObj *> &variables)
  : DensityObj(), _densities()
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(variables[i]->density());
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
        tmp.array() *= (-pdfs.col(j).array());
      } else {
        tmp.array() *= (1-cdfs.col(j).array());
      }
    }
    out += tmp;
  }
  out = -out;
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

void
MinimumDensityObj::sample(Eigen::VectorXd &out) const {
  Eigen::VectorXd tmp;
  Eigen::MatrixXd samples(out.size(), _densities.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->sample(tmp); samples.col(i) = tmp;
  }
  for (size_t i=0; i<out.size(); i++) {
    out[i] = samples.row(i).minCoeff();
  }
}


/* ********************************************************************************************* *
 * Implementation of MaximumObj
 * ********************************************************************************************* */
MaximumObj::MaximumObj(VarObj *a, VarObj *b)
  : VarObj(), _variables(), _density(0)
{
  if (! a->mutuallyIndep(b)) {
    AssumptionError err;
    err << "Cannot assemble maximum variable, arguments are not mutually independent.";
    throw err;
  }
  if (MaximumObj *max_a = dynamic_cast<MaximumObj *>(a)) {
    _variables = max_a->variables();
  } else {
    _variables.push_back(a);
  }
  if (MaximumObj *max_b = dynamic_cast<MaximumObj *>(b)) {
    for (size_t i=0; i<max_b->variables().size(); i++) {
      _variables.push_back(max_b->variables()[i]);
    }
  } else {
    _variables.push_back(b);
  }

  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    // Add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add variable itself
    _dependencies.insert(_variables[i]);
  }

  _density = new MaximumDensityObj(_variables);
}

MaximumObj::MaximumObj(const std::vector<VarObj *> &variables)
  : VarObj(), _variables(variables), _density(0)
{
  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    if (! this->mutuallyIndep(_variables[i])) {
      AssumptionError err;
      err << "Cannot assemble maximum variable, arguments are not mutually independent.";
      throw err;
    }
    // Add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add variable itself
    _dependencies.insert(_variables[i]);
  }

  _density = new MaximumDensityObj(_variables);
}

MaximumObj::~MaximumObj() {
  // pass...
}

void
MaximumObj::mark() {
  if (isMarked()) { return; }
  VarObj::mark();
  for (size_t i=0; i<_variables.size(); i++) {
    _variables[i]->mark();
  }
  _density->mark();
}

DensityObj *
MaximumObj::density() {
  return _density;
}


/* ********************************************************************************************* *
 * Implementation of MinimumObj
 * ********************************************************************************************* */
MinimumObj::MinimumObj(VarObj *a, VarObj *b)
  : VarObj(), _variables(), _density(0)
{
  if (! a->mutuallyIndep(b)) {
    AssumptionError err;
    err << "Cannot assemble minimum variable, arguments are not mutually independent.";
    throw err;
  }
  if (MinimumObj *max_a = dynamic_cast<MinimumObj *>(a)) {
    _variables = max_a->variables();
  } else {
    _variables.push_back(a);
  }
  if (MinimumObj *max_b = dynamic_cast<MinimumObj *>(b)) {
    for (size_t i=0; i<max_b->variables().size(); i++) {
      _variables.push_back(max_b->variables()[i]);
    }
  } else {
    _variables.push_back(b);
  }

  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    // Add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add variable itself
    _dependencies.insert(_variables[i]);
  }

  _density = new MinimumDensityObj(_variables);
}

MinimumObj::MinimumObj(const std::vector<VarObj *> &variables)
  : VarObj(), _variables(variables), _density(0)
{
  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    if (! this->mutuallyIndep(_variables[i])) {
      AssumptionError err;
      err << "Cannot assemble minimum variable, arguments are not mutually independent.";
      throw err;
    }
    // Add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add variable itself
    _dependencies.insert(_variables[i]);
  }

  _density = new MinimumDensityObj(_variables);
}

MinimumObj::~MinimumObj() {
  // pass...
}

void
MinimumObj::mark() {
  if (isMarked()) { return; }
  VarObj::mark();
  for (size_t i=0; i<_variables.size(); i++) {
    _variables[i]->mark();
  }
  _density->mark();
}

DensityObj *
MinimumObj::density() {
  return _density;
}

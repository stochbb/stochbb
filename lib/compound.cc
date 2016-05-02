#include "compound.hh"
#include "operators.hh"
#include "math.hh"
#include "logger.hh"
#include "distribution.hh"

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
 * Implementation of GenericCompoundObj
 * ********************************************************************************************* */
GenericCompoundObj::GenericCompoundObj(const std::vector<Var> &vars, const Distribution &distribution, const std::string &name)
  : CompoundObj(vars, name), _density(0)
{
  std::vector<DensityObj *> densities; densities.reserve(vars.size());
  for (size_t i=0; i<vars.size(); i++) {
    densities.push_back(*vars[i].density());
  }
  _density = new GenericCompoundDensityObj(*distribution, densities);
}

void
GenericCompoundObj::mark() {
  if (isMarked())
    return;
  CompoundObj::mark();
  if (_density)
    _density->mark();
}

Density
GenericCompoundObj::density() {
  return _density;
}

Distribution
GenericCompoundObj::distribution() {
  return _density->distribution();
}


/* ********************************************************************************************* *
 * Implementation of GenericCompoundDensityObj
 * ********************************************************************************************* */
GenericCompoundDensityObj::GenericCompoundDensityObj(DistributionObj *dist, const std::vector<DensityObj *> &params)
  : DensityObj(), _distribution(dist), _parameters(params)
{
  // pass...
}

GenericCompoundDensityObj::~GenericCompoundDensityObj() {
  // pass...
}

void
GenericCompoundDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  _distribution->mark();
  for (size_t i=0; i<_parameters.size(); i++) {
    _parameters[i]->mark();
  }
}

Distribution
GenericCompoundDensityObj::distribution() const {
  _distribution->ref();
  return _distribution;
}

size_t
GenericCompoundDensityObj::nParams() const {
  return _distribution->nParams();
}

Density
GenericCompoundDensityObj::parameter(size_t i) const {
  _parameters[i]->ref();
  return _parameters[i];
}

void
GenericCompoundDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  size_t Nstep = 100;
  double df = 1;

  // Get parameter values and PDFs
  Eigen::MatrixXd PDFs(Nstep, _parameters.size());
  Eigen::MatrixXd params(Nstep, _parameters.size());
  for (size_t i=0; i<_parameters.size(); i++) {
    double a, b;
    _parameters[i]->rangeEst(0.01, a, b);

    double x = a, dx = (b-a)/Nstep;
    df *= (b-a)/Nstep;

    for (size_t j=0; j<Nstep; j++, x+=dx)
      params(j,i) = x;

    _parameters[i]->eval(a,b, PDFs.col(i));
  }

  out.setZero();
  size_t N = std::pow(Nstep, _parameters.size());
  Eigen::VectorXd param(_parameters.size());
  std::vector<size_t> idxs(0, _parameters.size());
  double dx=(Tmax-Tmin)/out.size();

  // For each value in [Tmin, Tmax):
#pragma omp for
  for (int i=0; i<out.size(); i++) {
    double x = Tmin + i*dx;
    // "Integrate" over parameter space
    for (size_t j=0; j<N; j++) {
      // Get parameter indices
      _to_param_indices(j, Nstep, idxs);
      // Get parameter vector and prod of parameter PDFs
      double dp = 1;
      for (size_t k=0; k<_parameters.size(); k++) {
        param[k] = params(idxs[k], k);
        dp *= PDFs(idxs[k], k);
      }
      _distribution->pdf(x, x+dx, out.segment<1>(i), param);
      out(i) *= dp*df;
    }
  }
}

void
GenericCompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  size_t Nstep = 100;
  double df = 1;

  // Get parameter values and PDFs
  Eigen::MatrixXd PDFs(Nstep, _parameters.size());
  Eigen::MatrixXd params(Nstep, _parameters.size());
  for (size_t i=0; i<_parameters.size(); i++) {
    double a, b;
    _parameters[i]->rangeEst(0.01, a, b);

    double x = a, dx = (b-a)/Nstep;
    df *= (b-a)/Nstep;

    for (size_t j=0; j<Nstep; j++, x+=dx)
      params(j,i) = x;

    _parameters[i]->eval(a,b, PDFs.col(i));
  }

  out.setZero();
  size_t N = std::pow(Nstep, _parameters.size());
  Eigen::VectorXd param(_parameters.size());
  std::vector<size_t> idxs(0, _parameters.size());
  double dx=(Tmax-Tmin)/out.size();

  // For each value in [Tmin, Tmax):
#pragma omp for
  for (int i=0; i<out.size(); i++) {
    double x = Tmin + i*dx;
    // "Integrate" over parameter space
    for (size_t j=0; j<N; j++) {
      // Get parameter indices
      _to_param_indices(j, Nstep, idxs);
      // Get parameter vector and prod of parameter PDFs
      double dp = 1;
      for (size_t k=0; k<_parameters.size(); k++) {
        param[k] = params(idxs[k], k);
        dp *= PDFs(idxs[k], k);
      }
      _distribution->cdf(x, x+dx, out.segment<1>(i), param);
      out(i) *= dp*df;
    }
  }
}

void
GenericCompoundDensityObj::_to_param_indices(size_t i, size_t N, std::vector<size_t> &idxs) const {
  for (size_t j=0; j<idxs.size(); j++) {
    idxs[j] = i%N; i /= N;
  }
}

Density
GenericCompoundDensityObj::affine(double scale, double shift) const {
  GenericCompoundDensityObj *res = new GenericCompoundDensityObj(_distribution, _parameters);
  _distribution->affine(scale, shift, res->_parameters);
  return res;
}

void
GenericCompoundDensityObj::rangeEst(double alpha, double &a, double &b) const {
  Eigen::VectorXd A(_parameters.size()), B(_parameters.size());
  for (size_t i=0; i<_parameters.size(); i++) {
    _parameters[i]->rangeEst(alpha, A[i], B[i]);
  }
  double tmp;
  _distribution->quantile(a, tmp, alpha, A);
  _distribution->quantile(tmp, b, alpha, B);
}

int
GenericCompoundDensityObj::compare(const DensityObj &other) const {
  // Compare by density type
  if (int res = DensityObj::compare(other))
    return res;
  // same density class -> cast
  const GenericCompoundDensityObj &cother = static_cast<const GenericCompoundDensityObj &>(other);
  // Compare by distribution
  if (int res = _distribution->compare(*cother._distribution))
    return res;
  // Compare by parameter densities
  for (size_t i=0; i<_parameters.size(); i++) {
    if (int res = _parameters[i]->compare(*cother._parameters[i]))
      return res;
  }
  return 0;
}

void
GenericCompoundDensityObj::print(std::ostream &stream) const {
  stream << "<CompoundDensity of "; _distribution->print(stream);
  stream << " with [";
  _parameters[0]->print(stream);
  for (size_t i=1; i<_parameters.size(); i++) {
    stream << ", ";
    _parameters[i]->print(stream);
  }
  stream << "] #" << this << ">";
}




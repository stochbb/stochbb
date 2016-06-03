#include "compound.hh"
#include "operators.hh"
#include "math.hh"
#include "logger.hh"
#include "distribution.hh"
#include "reduction.hh"


using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of GenericCompoundObj
 * ********************************************************************************************* */
CompoundObj::CompoundObj(const std::vector<Var> &vars, const Distribution &distribution, const std::string &name)
  : DerivedVarObj(vars, name), _distribution(*distribution), _density(0), _parameters()
{
  std::vector<DensityObj *> densities;
  densities.reserve(vars.size()); _parameters.reserve(vars.size());
  for (size_t i=0; i<vars.size(); i++) {
    _parameters.push_back(*vars[i]);
    densities.push_back(*vars[i].density());
  }
  Density density = new CompoundDensityObj(distribution, densities);
  // try to simplify density
  while (CompoundReductionRule *rule = CompoundReductions::get().find(density)) {
    density = rule->apply(density);
  }
  // finally store density
  _density = *density;
}

void
CompoundObj::mark() {
  if (isMarked())
    return;
  DerivedVarObj::mark();
  if (_distribution)
    _distribution->mark();
  if (_density)
    _density->mark();
  for (size_t i=0; i<_parameters.size(); i++) {
    _parameters[i]->mark();
  }
}

Density
CompoundObj::density() {
  _density->ref();
  return _density;
}

Distribution
CompoundObj::distribution() {
  _distribution->ref();
  return _distribution;
}

Var
CompoundObj::parameter(size_t i) const {
  _parameters[i]->ref();
  return _parameters[i];
}

void
CompoundObj::print(std::ostream &stream) const {
  stream << "<Compound distr=";
  _distribution->print(stream);
  stream << " params=[";
  _parameters[0]->print(stream);
  for (size_t i=1; i<_parameters.size(); i++) {
    stream << ", "; _parameters[i]->print(stream);
  }
  stream << " desity="; _density->print(stream);
  stream << " #" << this << ">";
}

void
CompoundObj::sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                    Eigen::Ref<Eigen::MatrixXd> samples) const
{
  Eigen::VectorXd val(1);
  Eigen::VectorXd params(this->_parameters.size());
  for (int i=0; i<samples.rows(); i++) {
    // assemble parameter vector
    for (int j=0; j<params.size(); j++) {
      params(j) = samples(i, indices(j));
    }
    _distribution->sample(val, params);
    samples(i,outIdx) = val(0);
  }
}


/* ********************************************************************************************* *
 * Implementation of GenericCompoundDensityObj
 * ********************************************************************************************* */
CompoundDensityObj::CompoundDensityObj(const Distribution &dist, const std::vector<DensityObj *> &params)
  : DensityObj(), _distribution(*dist), _parameters(params)
{
  logDebug() << "Construct compound density object for distribution " << dist << "...";

  assume(_distribution->nParams() == _parameters.size());
}

CompoundDensityObj::~CompoundDensityObj() {
  // pass...
}

void
CompoundDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  if (_distribution)
    _distribution->mark();
  for (size_t i=0; i<_parameters.size(); i++) {
    _parameters[i]->mark();
  }
}

Distribution
CompoundDensityObj::distribution() const {
  // increment reference counter
  _distribution->ref();
  // transfer new reference to Distribution container
  return _distribution;
}

size_t
CompoundDensityObj::nParams() const {
  return _distribution->nParams();
}

Density
CompoundDensityObj::parameter(size_t i) const {
  // increment reference counter
  _parameters[i]->ref();
  // transfer new reference to Density container
  return _parameters[i];
}

void
CompoundDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  size_t Nstep = 100;
  double df = 1;

  // Get parameter values and PDFs of parameters
  Eigen::MatrixXd PDFs(Nstep, _parameters.size());
  Eigen::MatrixXd params(Nstep, _parameters.size());
  std::vector<size_t> Ns; Ns.reserve(_parameters.size());
  size_t N = 1;
  for (size_t i=0; i<_parameters.size(); i++) {
    AtomicDensityObj *p_atom = dynamic_cast<AtomicDensityObj *>(_parameters[i]);
    if (p_atom && dynamic_cast<DeltaDistributionObj *>(*p_atom->distribution())) {
      // Handle delta distributions
      Ns.push_back(1);
      params.col(i).setConstant(p_atom->parameter(0));
      PDFs.col(i).setConstant(1);
      N *= 1;
    } else {
      Ns.push_back(Nstep);
      N *= Nstep;
      double a, b;
      _parameters[i]->rangeEst(0.00001, a, b);
      double x = a, dx = (b-a)/Nstep; df *= dx;
      for (size_t j=0; j<Nstep; j++, x+=dx)
        params(j,i) = x;
      // get parameter PDF
      _parameters[i]->eval(a, b, PDFs.col(i));
    }
  }

  out.setZero();
  // the current parameter vector
  Eigen::VectorXd param(_parameters.size());
  // value indices for the i-th summand
  std::vector<size_t> idxs(_parameters.size());
  // temp vector holding the PDF for a specific param. vector
  Eigen::VectorXd tmp(out.size()); out.setZero();

  // For each value in [Tmin, Tmax):
  //   "Integrate" over parameter space
  for (size_t j=0; j<N; j++) {
    // Get parameter indices
    _to_param_indices(j, Ns, idxs);
    // Get parameter vector and prod of parameter PDFs
    double dp = 1;
    for (size_t k=0; k<_parameters.size(); k++) {
      param[k] = params(idxs[k], k);
      dp *= PDFs(idxs[k], k);
    }
    // eval distribution PDF for given parameter vector
    _distribution->pdf(Tmin, Tmax, tmp, param);
    // update result
    out += tmp*dp*df;
  }
}

void
CompoundDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  size_t Nstep = 100;
  double df = 1;

  // Get parameter values and PDFs of parameters
  Eigen::MatrixXd PDFs(Nstep, _parameters.size());
  Eigen::MatrixXd params(Nstep, _parameters.size());
  std::vector<size_t> Ns; Ns.reserve(_parameters.size());
  size_t N = 1;
  for (size_t i=0; i<_parameters.size(); i++) {
    AtomicDensityObj *p_atom = dynamic_cast<AtomicDensityObj *>(_parameters[i]);
    if (p_atom && dynamic_cast<DeltaDistributionObj *>(*p_atom->distribution())) {
      // Handle delta distributions
      Ns.push_back(1);
      params.col(i).setConstant(p_atom->parameter(0));
      PDFs.col(i).setConstant(1);
      N *= 1;
    } else {
      Ns.push_back(Nstep);
      N *= Nstep;
      double a, b;
      _parameters[i]->rangeEst(0.00001, a, b);
      double x = a, dx = (b-a)/Nstep; df *= dx;
      for (size_t j=0; j<Nstep; j++, x+=dx)
        params(j,i) = x;
      // get parameter PDF
      _parameters[i]->eval(a, b, PDFs.col(i));
    }
  }

  out.setZero();
  // the current parameter vector
  Eigen::VectorXd param(_parameters.size());
  // value indices for the i-th summand
  std::vector<size_t> idxs(_parameters.size());
  // temp vector holding the PDF for a specific param. vector
  Eigen::VectorXd tmp(out.size()); out.setZero();

  // For each value in [Tmin, Tmax):
  //   "Integrate" over parameter space
  for (size_t j=0; j<N; j++) {
    // Get parameter indices
    _to_param_indices(j, Ns, idxs);
    // Get parameter vector and prod of parameter PDFs
    double dp = 1;
    for (size_t k=0; k<_parameters.size(); k++) {
      param[k] = params(idxs[k], k);
      dp *= PDFs(idxs[k], k);
    }
    // eval distribution PDF for given parameter vector
    _distribution->cdf(Tmin, Tmax, tmp, param);
    // update result
    out += tmp*dp*df;
  }
}

void
CompoundDensityObj::_to_param_indices(size_t i, const std::vector<size_t> &Ns, std::vector<size_t> &idxs) const {
  for (size_t j=0; j<Ns.size(); j++) {
    idxs[j] = i%Ns[j]; i /= Ns[j];
  }
}

Density
CompoundDensityObj::affine(double scale, double shift) const {
  // increment reference counter
  _distribution->ref();
  // Copy, reference to _distribution is taken by the constructor
  CompoundDensityObj *res = new CompoundDensityObj(_distribution, _parameters);
  std::vector<Density> params; params.reserve(_parameters.size());
  for (size_t i=0; i<_parameters.size(); i++) {
    _parameters[i]->ref();
    params.push_back(Density(_parameters[i]));
  }
  // perform affine transform on parameter distributions
  _distribution->affine(scale, shift, params);
  // done, reference is transferred to Density container.
  return res;
}

void
CompoundDensityObj::rangeEst(double alpha, double &a, double &b) const {
  Eigen::VectorXd A(_parameters.size()), B(_parameters.size());
  for (size_t i=0; i<_parameters.size(); i++) {
    _parameters[i]->rangeEst(alpha, A[i], B[i]);
  }
  double tmp;
  _distribution->quantile(a, tmp, alpha, A);
  _distribution->quantile(tmp, b, alpha, B);
}

int
CompoundDensityObj::compare(const DensityObj &other) const {
  // Compare by density type
  if (int res = DensityObj::compare(other))
    return res;
  // same density class -> cast
  const CompoundDensityObj &cother = static_cast<const CompoundDensityObj &>(other);
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
CompoundDensityObj::print(std::ostream &stream) const {
  stream << "<CompoundDensity distr="; _distribution->print(stream);
  stream << " param=[";
  _parameters[0]->print(stream);
  for (size_t i=1; i<_parameters.size(); i++) {
    stream << ", ";
    _parameters[i]->print(stream);
  }
  stream << "] #" << this << ">";
}




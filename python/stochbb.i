%module stochbb

%{
#define SWIG_FILE_WITH_INIT
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "lib/api.hh"
#include "lib/randomvariable.hh"
#include "lib/density.hh"
%}

%include "numpy.i"

%init %{
import_array();
%}


%apply (double* INPLACE_ARRAY1, int DIM1) {(double* out, int N)}
%apply (double* INPLACE_FARRAY2, int DIM1, int DIM2) {(double* out, int Nrow, int Ncol)}

namespace sbb {

class Container
{
protected:
  Container();

  bool isNull() const;
};


class Density: public Container
{
protected:
  Density();
};


%extend Density {
  void eval(double Tmin, double Tmax, double* out, int N) const {
    Eigen::Map<Eigen::VectorXd> outMap(out, N);
    self->eval(Tmin, Tmax, outMap);
  }
}
%extend Density {
  void evalCDF(double Tmin, double Tmax, double* out, int N) const {
    Eigen::Map<Eigen::VectorXd> outMap(out, N);
    self->evalCDF(Tmin, Tmax, outMap);
  }
}


class Var: public Container
{
protected:
  Var();

public:
  Density density() const;
  bool dependsOn(const Var &var);
  bool mutuallyIndep(const Var &var);
};

%extend Var {
  Var __add__(Var *other) {
    return *self + *other;
  }
  Var __add__(double b) {
    return *self + b;
  }
  Var __radd__(double b) {
    return *self + b;
  }
  Var __mul__(double a) {
    return (*self)*a;
  }
  Var __rmul__(double a) {
    return (*self)*a;
  }
}

class DerivedVar: public Var
{
protected:
  DerivedVar();

public:
  size_t numVariables() const;
  Var variable(size_t i) const;
};


class ExactSampler: public Container
{
public:
  ExactSampler(const Var &X);
  ExactSampler(const Var &X1, const Var &X2);
  ExactSampler(const Var &X1, const Var &X2, const Var &X3);
  ExactSampler(const std::vector<Var> &variables);

  void sample(Eigen::MatrixXd &out) const;
};

%extend ExactSampler {
  void sample(double* out, int Nrow, int Ncol) {
    Eigen::Map<Eigen::MatrixXd> outMap(out, Nrow, Ncol);
    self->sample(outMap);
  }
}

class MarginalSampler: public Container
{
public:
  MarginalSampler(const Var &var, double Tmin, double Tmax, size_t steps);
};

%extend MarginalSampler {
  void sample(double* out, int N) {
    Eigen::Map<Eigen::VectorXd> outMap(out, N);
    self->sample(outMap);
  }
}


Var delta(double value);

Var uniform(double a, double b);

Var normal(double mu, double sigma);
Var normal(double mu, const Var &sigma);
Var normal(const Var &mu, double sigma);
Var normal(const Var &mu, const Var &sigma);

Var gamma(double k, double theta);
Var gamma(double k, const Var &theta);
Var gamma(const Var &k, double theta);
Var gamma(const Var &k, const Var &theta);

Var affine(const Var &var, double scale, double shift);

Var chain(const std::vector<Var> &vars);

Var minimum(const std::vector<Var> &variables);

Var maximum(const std::vector<Var> &variables);

%extend Container {
  bool isVar() const { return self->is<sbb::Var>(); }
  sbb::Var asVar() { return self->as<sbb::Var>(); }
  bool isDerivedVar() const { return self->is<sbb::DerivedVar>(); }
  sbb::DerivedVar asDerivedVar() { return self->as<sbb::DerivedVar>(); }
  bool isDensity() const { return self->is<sbb::Density>(); }
  sbb::Density asDensity() { return self->as<sbb::Density>(); }
}
}

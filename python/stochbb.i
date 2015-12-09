%module stochbb

%{
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "lib/api.hh"
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
import_array();
%}

namespace sbb {

class Container
{
protected:
  Container();
};


class Density: public Container
{
protected:
  Density();

public:
  void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
};


class Var: public Container
{
protected:
  Var();

public:
  Density density() const;
  bool dependsOn(const Var &var);
  bool mutuallyIndep(const Var &var);
};


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

}

%module stochbb

%{
#define SWIG_FILE_WITH_INIT
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "lib/stochbb.hh"
%}

%include "numpy.i"
%include "std_vector.i"
%include "std_string.i"
%include "typemaps.i"


%init %{
import_array();
%}


%apply (double* INPLACE_ARRAY1, int DIM1) {(double* out, int N)}
%apply (double* INPLACE_FARRAY2, int DIM1, int DIM2) {(double* out, int Nrow, int Ncol)}
%apply double *OUTPUT { double *aout, double *bout }


namespace stochbb {

class Error { };
%extend Error {
  std::string __str__() {
    return self->str();
  }
}

class Container
{
protected:
  Container();

  bool isNull() const;
};

%extend Container {
  std::string __str__() {
    std::stringstream buffer;
    buffer << *self;
    return buffer.str();
  }
}

class Density: public Container
{
protected:
  Density();
};


%extend Density {
  void eval(double Tmin, double Tmax, double* out, int N) const {
    Eigen::Map<Eigen::VectorXd> outMap(out, N);
    $self->eval(Tmin, Tmax, outMap);
  }

  void evalCDF(double Tmin, double Tmax, double* out, int N) const {
    Eigen::Map<Eigen::VectorXd> outMap(out, N);
    $self->evalCDF(Tmin, Tmax, outMap);
  }

  void rangeEst(double alpha, double *aout, double *bout) const {
    $self->rangeEst(alpha, *aout, *bout);
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

}

// Define vector-of-variables type:
//  This allows to pass list of variables to C++ functions taking
//  std::vector<Var> arguments.
namespace std {
  %template(varvector) vector<stochbb::Var>;
};

namespace stochbb {

class DerivedVar: public Var
{
protected:
  DerivedVar();

public:
  size_t numVariables() const;
  Var variable(size_t i) const;
};


class AffineTrafo: public DerivedVar
{
protected:
  AffineTrafo();

public:
  double scale() const;
  double shift() const;
};


class Chain: public DerivedVar
{
protected:
  Chain();

public:
  Chain(const std::vector<Var> &vars);
};


class Maximum: public DerivedVar
{
protected:
  Maximum();
};


class Minimum: public DerivedVar
{
protected:
  Minimum();
};


class Mixture: public DerivedVar
{
protected:
  Mixture();

public:
  double weight(size_t i) const;
};


class Conditional: public DerivedVar
{
protected:
  Conditional();
};


class CondChain: public DerivedVar
{
protected:
  CondChain();
};


class Compound: public DerivedVar
{
protected:
  Compound();
};


class ExactSampler: public Container
{
public:
  ExactSampler(const Var &X) throw (stochbb::TypeError);
  ExactSampler(const Var &X1, const Var &X2) throw (stochbb::TypeError);
  ExactSampler(const Var &X1, const Var &X2, const Var &X3) throw (stochbb::TypeError);
  ExactSampler(const std::vector<Var> &variables) throw (stochbb::TypeError);

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

Var weibull(double k, double lambda);
Var weibull(const Var& k, double lambda);
Var weibull(double k, const Var& lambda);
Var weibull(const Var& k, const Var& lambda);

Var affine(const Var &var, double scale, double shift);

Var chain(const std::vector<Var> &vars) throw( Error );

Var minimum(const Var &X1, const Var &X2) throw( Error );
Var minimum(const Var &X1, const Var &X2, const Var &X3) throw( Error );
Var minimum(const std::vector<Var> &variables) throw( Error );

Var maximum(const Var &X1, const Var &X2) throw( Error );
Var maximum(const Var &X1, const Var &X2, const Var &X3) throw( Error );
Var maximum(const std::vector<Var> &variables) throw( Error );

Var mixture(double w1, const Var &X1, double w2, const Var &X2) throw( Error );
Var mixture(double w1, const Var &X1, double w2, const Var &X2, double w3, const Var &X3) throw( Error );
Var mixture(const std::vector<double> weights, const std::vector<Var> &variables) throw( Error );

Var conditional(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2) throw( Error );
Var condchain(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2) throw( Error );

double dnorm(double x);
double pnorm(double x);
double qnorm(double p);

double dgamma(double x, double k, double theta);
double pgamma(double x, double k, double theta);
double qgamma(double p, double k, double theta);

double dweibull(double x, double k, double theta);
double pweibull(double x, double k, double theta);
double qweibull(double p, double k, double theta);

%extend Container {
  bool isDensity() const { return self->is<stochbb::Density>(); }
  stochbb::Density asDensity() { return self->as<stochbb::Density>(); }

  bool isVar() const { return self->is<stochbb::Var>(); }
  stochbb::Var asVar() { return self->as<stochbb::Var>(); }

  bool isDerivedVar() const { return self->is<stochbb::DerivedVar>(); }
  stochbb::DerivedVar asDerivedVar() { return self->as<stochbb::DerivedVar>(); }

  bool isAffineTrafo() const { return self->is<stochbb::AffineTrafo>(); }
  stochbb::AffineTrafo asAffineTrafo() { return self->as<stochbb::AffineTrafo>(); }

  bool isChain() const { return self->is<stochbb::Chain>(); }
  stochbb::Chain asChain() { return self->as<stochbb::Chain>(); }

  bool isMaximum() const { return self->is<stochbb::Maximum>(); }
  stochbb::Maximum asMaximum() { return self->as<stochbb::Maximum>(); }

  bool isMinimum() const { return self->is<stochbb::Minimum>(); }
  stochbb::Minimum asMinimum() { return self->as<stochbb::Minimum>(); }

  bool isMixture() const { return self->is<stochbb::Mixture>(); }
  stochbb::Mixture asMixture() { return self->as<stochbb::Mixture>(); }

  bool isConditional() const { return self->is<stochbb::Conditional>(); }
  stochbb::Conditional asConditional() { return self->as<stochbb::Conditional>(); }

  bool isCompound() const { return self->is<stochbb::Compound>(); }
  stochbb::Compound asCompound() { return self->as<stochbb::Compound>(); }
}

class LogMessage
{
public:
  typedef enum {
    DEBUG = 0,
    INFO,
    WARNING,
    ERROR
  } Level;

public:
  LogMessage();
  LogMessage(const std::string &filename, int line, Level level, const std::string &message);
  LogMessage(const LogMessage &other);

  const std::string &filename() const;
  int linenumber() const;
  Level level() const;
  const std::string &message() const;
  const std::time_t &timestamp() const;
};


class LogHandler: public Container
{
protected:
  LogHandler(LogHandlerObj *obj);

public:
  void handleMessage(const LogMessage &msg);
};


class IOLogHandler: public LogHandler
{
public:
  IOLogHandler(std::ostream &stream=std::cerr, LogMessage::Level level=LogMessage::DEBUG);
};

class Logger
{
protected:
  Logger();

public:
  /** Logs a message. */
  static void log(const LogMessage &msg);
  /** Adds a handler to the logger, the ownership is transferred to the @c Logger. */
  static void addHandler(const LogHandler &handler);
};

}

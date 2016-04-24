#include <RcppCommon.h>

#include <stochbb/stochbb.hh>
using namespace stochbb;

RCPP_EXPOSED_CLASS_NODECL(Container);
RCPP_EXPOSED_CLASS_NODECL(Density);
RCPP_EXPOSED_CLASS_NODECL(Var);

#include <RcppEigen.h>

using namespace Rcpp;

typedef Eigen::Map<Eigen::VectorXd> MapVec;


std::string printContainer(Container *cont) {
  std::stringstream stream;
  (*cont)->print(stream);
  return stream.str();
}


RCPP_MODULE(stochbb) {
  class_<Container>("Container")
      .property("isNull", &Container::isNull)
      .method("repr", &printContainer);

  class_<Density>("Density")
      .derives<Container>("Container")
      .method<void, double, double, MapVec>("eval", (void (Density::*)(double, double, MapVec)) &Density::eval)
      .method<void, double, double, MapVec>("evalCDF", (void (Density::*)(double, double, MapVec)) &Density::evalCDF);

  class_<Var>("Var")
      .derives<Container>("Container")
      .property("density", &Var::density);

  function<Var, double>("delta", &delta);
  function<Var, double, double, const std::string &>(
        "uniform", &uniform,
        List::create( _["a"], _["b"], _["name"]=""));

  function<Var, double, double, const std::string &>(
        "normal", &normal,
        List::create( _["mu"], _["sigma"], _["name"]=""));
  function<Var, double, const Var &, const std::string &>(
        "normal", &normal,
        List::create( _["mu"], _["sigma"], _["name"]=""));
  function<Var, const Var &, double, const std::string &>(
        "normal", &normal,
        List::create( _["mu"], _["sigma"], _["name"]=""));
  function<Var, const Var &, const Var &, const std::string &>(
        "normal", &normal,
        List::create( _["mu"], _["sigma"], _["name"]=""));

  function<Var, double, double, const std::string &>(
        "gamma", &gamma,
        List::create( _["k"], _["theta"], _["name"]=""));
  function<Var, double, const Var &, const std::string &>(
        "gamma", &gamma,
        List::create( _["k"], _["theta"], _["name"]=""));
  function<Var, const Var &, double, const std::string &>(
        "gamma", &gamma,
        List::create( _["k"], _["theta"], _["name"]=""));
  function<Var, const Var &, const Var &, const std::string &>(
        "gamma", &gamma,
        List::create( _["k"], _["theta"], _["name"]=""));

  function<Var, double, double, const std::string &>(
        "invgamma", &invgamma,
        List::create( _["alpha"], _["beta"], _["name"]=""));
  function<Var, double, const Var &, const std::string &>(
        "invgamma", &invgamma,
        List::create( _["alpha"], _["beta"], _["name"]=""));
  function<Var, const Var &, double, const std::string &>(
        "invgamma", &invgamma,
        List::create( _["alpha"], _["beta"], _["name"]=""));
  function<Var, const Var &, const Var &, const std::string &>(
        "invgamma", &invgamma,
        List::create( _["alpha"], _["beta"], _["name"]=""));

  function<Var, double, double, const std::string &>(
        "weibull", &weibull,
        List::create( _["k"], _["lambda"], _["name"]=""));
  function<Var, double, const Var &, const std::string &>(
        "weibull", &weibull,
        List::create( _["k"], _["lambda"], _["name"]=""));
  function<Var, const Var &, double, const std::string &>(
        "weibull", &weibull,
        List::create( _["k"], _["lambda"], _["name"]=""));
  function<Var, const Var &, const Var &, const std::string &>(
        "weibull", &weibull,
        List::create( _["k"], _["lambda"], _["name"]=""));

  function<Var, const Var &, double, double>(
        "affine", &affine);

  function<Var, const Var &, const Var &>(
        "chain", (Var (*)(const Var &, const Var &)) &chain);
  function<Var, const Var &, const Var &, const Var &>(
        "chain", (Var (*)(const Var &, const Var &, const Var &)) &chain);

  function<Var, const Var &, const Var &>(
        "minimum", (Var (*)(const Var &, const Var &)) &minimum);
  function<Var, const Var &, const Var &, const Var &>(
        "minimum", (Var (*)(const Var &, const Var &, const Var &)) &minimum);

  function<Var, const Var &, const Var &>(
        "maximum", (Var (*)(const Var &, const Var &)) &maximum);
  function<Var, const Var &, const Var &, const Var &>(
        "maximum", (Var (*)(const Var &, const Var &, const Var &)) &maximum);

  function<Var, double, const Var &, double, const Var &>(
        "mixture", (Var (*)(double, const Var &, double, const Var &)) &mixture);
  function<Var, double, const Var &, double, const Var &, double, const Var &>(
        "mixture", (Var (*)(double, const Var &, double, const Var &, double, const Var &)) &mixture);

  Var mixture(double wX1, const Var &X1, double wX2, const Var &X2, double wX3, const Var &X3);

  function<Var, const Var &, const Var &, const Var &, const Var &>(
        "conditional", &conditional);

  function<Var, const Var &, const Var &, const Var &, const Var &>(
        "condchain", &stochbb::condchain);
}

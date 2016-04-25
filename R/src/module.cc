#include <RcppCommon.h>

#include <stochbb/stochbb.hh>
using namespace stochbb;

RCPP_EXPOSED_CLASS_NODECL(Container);
RCPP_EXPOSED_CLASS_NODECL(Density);
RCPP_EXPOSED_CLASS_NODECL(Var);

#include <RcppEigen.h>

using namespace Rcpp;

typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;


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
      .property("density", &Var::density)
      .property("name", &Var::name, &Var::setName)
      .method("mutuallyIndep", &Var::mutuallyIndep)
      .method("dependsOn", &Var::dependsOn);

  class_<ExactSampler>("ExactSampler")
      .derives<Container>("Container")
      .constructor<const Var &>()
      .constructor<const Var &, const Var &>()
      .constructor<const Var &, const Var &, const Var &>()
      .method<void, MapMat>("sample", (void (ExactSampler::*)(MapMat)) &ExactSampler::sample);

  function<bool, const Var &, const Var &>(
        "independent", (bool (*)(const Var &, const Var &)) &independent);
  function<bool, const Var &, const Var &, const Var &>(
        "independent", (bool (*)(const Var &, const Var &, const Var &)) &independent);

  function<Var, double>("delta", &delta);

  function<Var, double, double, const std::string &>(
        "uniform", &uniform);

  function<Var, double, double, const std::string &>(
        "normalrv", &normal);
  function<Var, const Var &, const Var &, const std::string &>(
        "compnormalrv", &normal);

  function<Var, double, double, const std::string &>(
        "gammarv", &gamma);
  function<Var, const Var &, const Var &, const std::string &>(
        "compgammarv", &gamma);

  function<Var, double, double, const std::string &>(
        "invgammarv", &invgamma);
  function<Var, const Var &, const Var &, const std::string &>(
        "compinvgammarv", &invgamma);

  function<Var, double, double, const std::string &>(
        "weibullrv", &weibull);
  function<Var, const Var &, const Var &, const std::string &>(
        "compweibullrv", &weibull);

  function<Var, const Var &, double, double>(
        "affine", &affine);

  function<Var, const Var &, const Var &>(
        "chain", (Var (*)(const Var &, const Var &)) &chain);

  function<Var, const Var &, const Var &>(
        "minimum", (Var (*)(const Var &, const Var &)) &minimum);

  function<Var, const Var &, const Var &>(
        "maximum", (Var (*)(const Var &, const Var &)) &maximum);

  function<Var, double, const Var &, double, const Var &>(
        "mixture", (Var (*)(double, const Var &, double, const Var &)) &mixture);

  function<Var, const Var &, const Var &, const Var &, const Var &>(
        "conditional", &conditional);

  function<Var, const Var &, const Var &, const Var &, const Var &>(
        "condchain", &stochbb::condchain);
}

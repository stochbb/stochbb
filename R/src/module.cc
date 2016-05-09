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

std::string print_Container(Container *cont) {
  std::stringstream stream;
  (*cont)->print(stream);
  return stream.str();
}

Rcpp::List Density_rangeEst(Density *density, double alpha) {
  double a, b;
  density->rangeEst(alpha, a, b);
  Rcpp::List res; res.push_back(a); res.push_back(b);
  return res;
}

ExactSampler* new_ExactSampler(const Rcpp::List &args) {
  std::vector<Var> vars;
  vars.reserve(args.size());
  for (Rcpp::List::const_iterator item=args.begin(); item != args.end(); item++) {
    vars.push_back(Rcpp::as<Var>(*item));
  }
  return new ExactSampler(vars);
}

void ExactSamper_sample(ExactSampler *sampler, Eigen::Map<Eigen::MatrixXd> out) {
  sampler->sample(out);
}

bool independent_List(const Rcpp::List &args) {
  std::vector<Var> vars;
  vars.reserve(args.size());
  for (Rcpp::List::const_iterator item=args.begin(); item != args.end(); item++) {
    vars.push_back(Rcpp::as<Var>(*item));
  }
  return independent(vars);
}

Var chain_List(const Rcpp::List &args) {
  std::vector<Var> vars;
  vars.reserve(args.size());
  for (Rcpp::List::const_iterator item=args.begin(); item != args.end(); item++) {
    vars.push_back(Rcpp::as<Var>(*item));
  }
  return chain(vars);
}

Var minimum_List(const Rcpp::List &args) {
  std::vector<Var> vars;
  vars.reserve(args.size());
  for (Rcpp::List::const_iterator item=args.begin(); item != args.end(); item++) {
    vars.push_back(Rcpp::as<Var>(*item));
  }
  return minimum(vars);
}

Var maximum_List(const Rcpp::List &args) {
  std::vector<Var> vars;
  vars.reserve(args.size());
  for (Rcpp::List::const_iterator item=args.begin(); item != args.end(); item++) {
    vars.push_back(Rcpp::as<Var>(*item));
  }
  return maximum(vars);
}


RCPP_MODULE(stochbb) {
  class_<Container>("Container")
      .property("isNull", &Container::isNull)
      .method("repr", &print_Container);

  class_<Density>("Density")
      .derives<Container>("Container")
      .method<void, double, double, MapVec>("eval", (void (Density::*)(double, double, MapVec)) &Density::eval)
      .method<void, double, double, MapVec>("evalCDF", (void (Density::*)(double, double, MapVec)) &Density::evalCDF)
      .method<Rcpp::List, double>("rangeEst", &Density_rangeEst);

  class_<Var>("Var")
      .derives<Container>("Container")
      .property("density", &Var::density)
      .property("name", &Var::name, &Var::setName)
      .method("mutuallyIndep", &Var::mutuallyIndep)
      .method("dependsOn", &Var::dependsOn);

  class_<ExactSampler>("ExactSampler")
      .derives<Container>("Container")
      .factory<const List &>(&new_ExactSampler)
      .method<void, MapMat>("sample", &ExactSamper_sample);

  function<bool, const List &>(
        "_independent", &independent_List);

  function<Var, double>("delta", &delta);

  function<Var, double, double, const std::string &>(
        "_uniformrv", &uniform);

  function<Var, double, double, const std::string &>(
        "_normalrv", &normal);
  function<Var, const Var &, const Var &, const std::string &>(
        "_compnormalrv", &normal);

  function<Var, double, double, const std::string &>(
        "_gammarv", &gamma);
  function<Var, const Var &, const Var &, const std::string &>(
        "_compgammarv", &gamma);

  function<Var, double, double, const std::string &>(
        "_invgammarv", &invgamma);
  function<Var, const Var &, const Var &, const std::string &>(
        "_compinvgammarv", &invgamma);

  function<Var, double, double, const std::string &>(
        "_weibullrv", &weibull);
  function<Var, const Var &, const Var &, const std::string &>(
        "_compweibullrv", &weibull);

  function<Var, double, const std::string &>(
        "_studtrv", &studt);
  function<Var, const Var &, const std::string &>(
        "_compstudtrv", &studt);

  function<Var, const Var &, double, double>(
        "affine", &affine);

  function<Var, const List &>(
        "_chain", &chain_List);

  function<Var, const List &>(
        "_minimum", &minimum_List);

  function<Var, const List &>(
        "_maximum", &maximum_List);

  function<Var, double, const Var &, double, const Var &>(
        "mixture", (Var (*)(double, const Var &, double, const Var &)) &mixture);

  function<Var, const Var &, const Var &, const Var &, const Var &>(
        "conditional", &conditional);

  function<Var, const Var &, const Var &, const Var &, const Var &>(
        "condchain", &stochbb::condchain);
}

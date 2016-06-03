#include "operators.hh"
#include "randomvariable.hh"
#include "affinetrafo.hh"
#include "chain.hh"
#include "minmax.hh"
#include "mixture.hh"
#include "conditional.hh"
#include "compound.hh"
#include "distribution.hh"
#include "logger.hh"

using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of sbb::delta()
 * ********************************************************************************************* */
bool
stochbb::is_delta(const Density &dens) {
  AtomicDensity at = dens.as<AtomicDensity>();
  return at.isValid() && at.distribution().is<DeltaDistribution>();
}

Var
stochbb::delta(double value, const std::__1::string &name) {
  Eigen::VectorXd param(1); param << value;
  return new AtomicVarObj(new AtomicDensityObj(new DeltaDistributionObj(), param), name);
}

/* ********************************************************************************************* *
 * Implementation of sbb::uniform()
 * ********************************************************************************************* */
bool
stochbb::is_uniform(const Density &dens) {
  AtomicDensity at = dens.as<AtomicDensity>();
  return at.isValid() && at.distribution().is<UniformDistribution>();
}

Var
stochbb::uniform(double a, double b, const std::string &name) throw (Error) {
  Eigen::VectorXd param(2); param << a, b;
  return new AtomicVarObj(new AtomicDensityObj(new UniformDistributionObj(), param), name);
}

/* ********************************************************************************************* *
 * Implementation of sbb::normal()
 * ********************************************************************************************* */
bool
stochbb::is_normal(const Density &dens) {
  AtomicDensity at = dens.as<AtomicDensity>();
  return at.isValid() && at.distribution().is<NormalDistribution>();
}

Var
stochbb::normal(double mu, double sigma, const std::string &name) throw (Error) {
  Eigen::VectorXd param(2); param << mu, sigma;
  return new AtomicVarObj(new AtomicDensityObj(new NormalDistributionObj(), param), name);
}

Var
stochbb::normal(const Var &mu, double sigma, const std::string &name) throw (Error) {
  // If mu is delta distributed -> simplify to atomic random variable
  AtomicDensity mu_atomic = mu.density().as<AtomicDensity>();
  if (mu_atomic.isValid() && mu_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::normal(mu_atomic.parameter(0), sigma, name);
  }
  // Otherwise assemble CompoundVar
  std::vector<Var> vars = {mu, delta(sigma)};
  return new CompoundObj(vars, new NormalDistributionObj(), name);
}

Var
stochbb::normal(double mu, const Var &sigma, const std::string &name) throw (Error) {
  // If sigma is delta distributed -> simplify to atomic random variable
  AtomicDensity sigma_atomic = sigma.density().as<AtomicDensity>();
  if (sigma_atomic.isValid() && dynamic_cast<DeltaDistributionObj *>(*sigma_atomic->distribution())) {
    return stochbb::normal(mu, sigma_atomic.parameter(0), name);
  }
  // Otherwise assemble CompoundVar
  std::vector<Var> vars = {delta(mu), sigma};
  return new CompoundObj(vars, new NormalDistributionObj(), name);
}

Var
stochbb::normal(const Var &mu, const Var &sigma, const std::string &name) throw (Error) {
  AtomicDensity mu_atomic = mu.density().as<AtomicDensity>();
  AtomicDensity sigma_atomic = sigma.density().as<AtomicDensity>();
  // If mu is delta distributed -> simplify to atomic random variable
  if (mu_atomic.isValid() && mu_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::normal(mu_atomic.parameter(0), sigma, name);
  }
  // else if sigma is delta distributed -> simplify to atomic random variable
  if (sigma_atomic.isValid() && sigma_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::normal(mu, sigma_atomic.parameter(0), name);
  }
  // Otherwise assemble CompoundVar
  std::vector<Var> vars = {mu, sigma};
  return new CompoundObj(vars, new NormalDistributionObj(), name);
}


/* ********************************************************************************************* *
 * Implementation of sbb::gamma()
 * ********************************************************************************************* */
bool
stochbb::is_gamma(const Density &dens) {
  AtomicDensity at = dens.as<AtomicDensity>();
  return at.isValid() && at.distribution().is<GammaDistribution>();
}

Var
stochbb::gamma(double k, double theta, const std::string &name) throw (Error) {
  Eigen::VectorXd param(3); param << k, theta, 0;
  return new AtomicVarObj(new AtomicDensityObj(new GammaDistributionObj(), param), name);
}

Var
stochbb::gamma(const Var &k, double theta, const std::string &name) throw (Error) {
  // If mu is delta distributed -> simplify to atomic random variable
  AtomicDensity k_atomic = k.density().as<AtomicDensity>();
  if (k_atomic.isValid() && k_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::gamma(k_atomic.parameter(0), theta, name);
  }
  // Otherwise assemble CompoundVar
  std::vector<Var> vars = {k, delta(theta), delta(0)};
  return new CompoundObj(vars, new GammaDistributionObj(), name);
}

Var
stochbb::gamma(double k, const Var &theta, const std::string &name) throw (Error) {
  // If sigma is delta distributed -> simplify to atomic random variable
  AtomicDensity theta_atomic = theta.density().as<AtomicDensity>();
  if (theta_atomic.isValid() && theta_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::gamma(k, theta_atomic.parameter(0), name);
  }
  // Otherwise assemble CompoundVar
  std::vector<Var> vars = {delta(k), theta, delta(0)};
  return new CompoundObj(vars, new GammaDistributionObj(), name);
}

Var
stochbb::gamma(const Var &k, const Var &theta, const std::string &name) throw (Error) {
  AtomicDensity k_atomic = k.density().as<AtomicDensity>();
  AtomicDensity theta_atomic = theta.density().as<AtomicDensity>();
  // If k is delta distributed -> simplify to atomic random variable
  if (k_atomic.isValid() && k_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::gamma(k_atomic.parameter(0), theta, name);
  }
  // else if theta is delta distributed -> simplify to atomic random variable
  if (theta_atomic.isValid() && theta_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::gamma(k, theta_atomic.parameter(0), name);
  }
  // Otherwise assemble CompoundVar
  std::vector<Var> vars = {k, theta, delta(0)};
  return new CompoundObj(vars, new GammaDistributionObj(), name);
}


/* ********************************************************************************************* *
 * Implementation of sbb::invgamma()
 * ********************************************************************************************* */
bool
stochbb::is_invgamma(const Density &dens) {
  AtomicDensity at = dens.as<AtomicDensity>();
  return at.isValid() && at.distribution().is<InvGammaDistribution>();
}

Var
stochbb::invgamma(double alpha, double beta, const std::string &name) throw (Error) {
  Eigen::VectorXd param(3); param << alpha, beta, 0;
  return new AtomicVarObj(new AtomicDensityObj(new InvGammaDistributionObj(), param), name);
}

Var
stochbb::invgamma(const Var &alpha, double beta, const std::string &name) throw (Error) {
  // If alpha is delta distributed -> simplify to atomic random variable
  AtomicDensity alpha_atomic = alpha.density().as<AtomicDensity>();
  if (alpha_atomic.isValid() && alpha_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::invgamma(alpha_atomic.parameter(0), beta, name);
  }
  // Otherwise assemble CompoundVar
  return new CompoundObj(std::vector<Var> {alpha, delta(beta), delta(0)},
                         new InvGammaDistributionObj(), name);
}

Var
stochbb::invgamma(double alpha, const Var &beta, const std::string &name) throw (Error) {
  // If beta is delta distributed -> simplify to atomic random variable
  AtomicDensity beta_atomic = beta.density().as<AtomicDensity>();
  if (beta_atomic.isValid() && beta_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::invgamma(alpha, beta_atomic.parameter(0), name);
  }
  // Otherwise assemble CompoundVar
  return new CompoundObj(std::vector<Var> {delta(alpha), beta, delta(0)},
                         new InvGammaDistributionObj(), name);
}

Var
stochbb::invgamma(const Var &alpha, const Var &beta, const std::string &name) throw (Error) {
  AtomicDensity alpha_atomic = alpha.density().as<AtomicDensity>();
  AtomicDensity beta_atomic = beta.density().as<AtomicDensity>();
  // If alpha is delta distributed -> simplify to atomic random variable
  if (alpha_atomic.isValid() && alpha_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::invgamma(alpha_atomic.parameter(0), beta, name);
  }
  // else if beta is delta distributed -> simplify to atomic random variable
  if (beta_atomic.isValid() && beta_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::invgamma(alpha, beta_atomic.parameter(0), name);
  }
  // Otherwise assemble CompoundVar
  return new CompoundObj(std::vector<Var> {alpha, beta, delta(0)},
                         new InvGammaDistributionObj(), name);
}


/* ********************************************************************************************* *
 * Implementation of sbb::weibull()
 * ********************************************************************************************* */
bool
stochbb::is_weibull(const Density &dens) {
  AtomicDensity at = dens.as<AtomicDensity>();
  return at.isValid() && at.distribution().is<WeibullDistribution>();
}

Var
stochbb::weibull(double k, double lambda, const std::string &name) throw (Error) {
  Eigen::VectorXd param(3); param << k, lambda, 0;
  return new AtomicVarObj(new AtomicDensityObj(new WeibullDistributionObj(), param), name);
}

Var
stochbb::weibull(const Var &k, double lambda, const std::string &name) throw (Error) {
  AtomicDensity k_atomic = k.density().as<AtomicDensity>();
  if (k_atomic.isValid() && k_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::weibull(k_atomic.parameter(0), lambda, name);
  }
  return new CompoundObj(std::vector<Var> {k, delta(lambda), delta(0)},
                         new WeibullDistributionObj(), name);
}

Var
stochbb::weibull(double k, const Var &lambda, const std::string &name) throw (Error) {
  AtomicDensity lambda_atomic = lambda.density().as<AtomicDensity>();
  if (lambda_atomic.isValid() && lambda_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::weibull(k, lambda_atomic.parameter(0), name);
  }
  return new CompoundObj(std::vector<Var> {delta(k), lambda, delta(0)},
                         new WeibullDistributionObj(), name);
}

Var
stochbb::weibull(const Var &k, const Var& lambda, const std::string &name) throw (Error) {
  AtomicDensity k_atomic = k.density().as<AtomicDensity>();
  AtomicDensity lambda_atomic = lambda.density().as<AtomicDensity>();
  if (k_atomic.isValid() && k_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::weibull(k_atomic.parameter(0), lambda, name);
  }
  if (lambda_atomic.isValid() && lambda_atomic.distribution().is<DeltaDistribution>()) {
    return stochbb::weibull(k, lambda_atomic.parameter(0), name);
  }
  return new CompoundObj(std::vector<Var> {k, lambda, delta(0)},
                         new WeibullDistributionObj(), name);
}

/* ********************************************************************************************* *
 * Implementation of sbb::studt()
 * ********************************************************************************************* */
bool
stochbb::is_studt(const Density &dens) {
  AtomicDensity at = dens.as<AtomicDensity>();
  return at.isValid() && at.distribution().is<StudtDistribution>();
}

Var
stochbb::studt(double nu, const std::string &name) throw (Error) {
  Eigen::VectorXd param(3); param << nu, 1, 0;
  return new AtomicVarObj(new AtomicDensityObj(new StudtDistributionObj(), param), name);
}

Var
stochbb::studt(const Var &nu, const std::string &name) throw (Error) {
  AtomicDensity nu_a = nu.density().as<AtomicDensity>();
  if (nu_a.isValid() && nu_a.distribution().is<DeltaDistribution>()) {
    return studt(nu_a.parameter(0), name);
  }
  return new CompoundObj(std::vector<Var> {nu, delta(1), delta(0)},
                         new StudtDistributionObj(), name);
}


/* ********************************************************************************************* *
 * Implementation of minimum
 * ********************************************************************************************* */
Var stochbb::minimum(const Var &X1, const Var &X2) throw (Error) {
  return minimum(std::vector<Var> {X1, X2});
}

Var stochbb::minimum(const Var &X1, const Var &X2, const Var &X3) throw (Error) {
  return minimum(std::vector<Var> {X1, X2, X3});
}

Var
stochbb::minimum(const std::vector<Var> &variables) throw (Error) {
  std::vector<Var> args(variables);
  Var common = splitCommon(args);
  if (! independent(args)) {
    AssumptionError err;
    err << "Cannot construct maximum of variables, variables are not independent.";
    throw err;
  }
  return Var(new MinimumObj(args)) + common;
}


/* ********************************************************************************************* *
 * Implementation of maximum
 * ********************************************************************************************* */
Var stochbb::maximum(const Var &X1, const Var &X2) throw (Error) {
  return maximum(std::vector<Var> {X1, X2});
}

Var stochbb::maximum(const Var &X1, const Var &X2, const Var &X3) throw (Error) {
  return maximum(std::vector<Var> {X1, X2, X3});
}

Var
stochbb::maximum(const std::vector<Var> &variables) throw (Error) {
  std::vector<Var> args(variables);
  Var common = splitCommon(args);
  if (! independent(args)) {
    AssumptionError err;
    err << "Cannot construct maximum of variables, variables are not independent.";
    throw err;
  }
  return Var(new MaximumObj(args)) + common;
}


/* ********************************************************************************************* *
 * Implementation of independent
 * ********************************************************************************************* */
bool
stochbb::independent(const std::vector<Var> &vars) {
  for (size_t i=0; i<vars.size(); i++) {
    for (size_t j=(i+1); j<vars.size(); j++) {
      if (! vars[i].mutuallyIndep(vars[j])) {
        return false;
      }
    }
  }
  return true;
}

bool
stochbb::independent(const Var &a, const Var &b) {
  return a.mutuallyIndep(b);
}

bool
stochbb::independent(const Var &a, const Var &b, const Var &c) {
  return independent(std::vector<Var> {a, b, c});
}

/* ********************************************************************************************* *
 * Implementation of splitCommon
 * ********************************************************************************************* */
Var
stochbb::splitCommon(std::vector<Var> &vars) {
  if (0 == vars.size())
    return delta(0);
  if (1 == vars.size()) {
    Var res = vars.back();
    vars.pop_back();
    return res;
  }

  // Find the common part of the variables formed as a sum
  std::vector<VarSet> indepvars; indepvars.reserve(vars.size());
  for (size_t i=0; i<vars.size(); i++) {
    if (vars[i].is<Chain>()) {
      indepvars.push_back(VarSet());
      Chain chain = vars[i].as<Chain>();
      for (size_t i=0; i<chain.numVariables(); i++) {
        indepvars.back().add(chain.variable(i));
      }
    } else {
      indepvars.push_back(VarSet());
      indepvars.back().add(vars[i]);
    }
  }

  // Get the intersection of all variable sets
  VarSet common = indepvars[0].intersect(indepvars[1]);
  for (size_t j=2; j<indepvars.size(); j++) {
    common = common.intersect(indepvars[j]);
  }

  // Remove common part from all sets
  for (size_t i=0; i<indepvars.size(); i++) {
    indepvars[i] = indepvars[i].difference(common);
  }

  // Assemble result
  for (size_t i=0; i<indepvars.size(); i++) {
    if (1 == indepvars[i].size()) {
      vars[i] = *indepvars[i].begin();
    } else {
      std::vector<Var> chain_args; chain_args.reserve(indepvars[i].size());
      VarSet::iterator item = indepvars[i].begin();
      for (; item != indepvars[i].end(); item++) {
        chain_args.push_back(*item);
      }
      vars[i] = new ChainObj(chain_args);
    }
  }

  if (common.isEmpty())
    return delta(0);
  if (1 == common.size())
    return *common.begin();

  std::vector<Var> chain_args; chain_args.reserve(common.size());
  VarSet::iterator item = common.begin();
  for (; item != common.end(); item++) {
    chain_args.push_back(*item);
  }

  return new ChainObj(chain_args);
}


/* ********************************************************************************************* *
 * Implementation of chain
 * ********************************************************************************************* */
Var
stochbb::chain(const std::vector<Var> &vars) throw (Error) {
  std::vector<Var> variables;
  variables.reserve(2*vars.size());

  // Flatten chains
  for (size_t i=0; i<vars.size(); i++) {
    if (vars[i].is<Chain>()) {
      Chain chain = vars[i].as<Chain>();
      for (size_t j=0; j<chain.numVariables(); j++) {
        variables.push_back(chain.variable(j));
      }
    } else {
      variables.push_back(vars[i]);
    }
  }

  // Check for mutual independence
  if (! independent(variables)) {
    AssumptionError err;
    err << "Cannot construct chain from mutually depending random variables.";
    throw err;
  }

  return Chain(variables);
}

Var
stochbb::chain(const Var &X1, const Var &X2) throw (Error) {
  return chain(std::vector<Var> {X1, X2});
}

Var
stochbb::chain(const Var &X1, const Var &X2, const Var &X3) throw (Error) {
  return chain(std::vector<Var> {X1, X2, X3});
}


/* ********************************************************************************************* *
 * Implementation of affine
 * ********************************************************************************************* */
Var
stochbb::affine(const Var &var, double scale, double shift) {
  // Flatten affine trafo objects
  if (var.is<AffineTrafo>()) {
    AffineTrafo a = var.as<AffineTrafo>();
    return new AffineTrafoObj(scale*a.scale(), scale*a.shift()+shift, a.variable(0));
  }
  return new AffineTrafoObj(scale, shift, var);
}


/* ********************************************************************************************* *
 * Implementation of mixture
 * ********************************************************************************************* */
Var
stochbb::mixture(double wX1, const Var &X1, double wX2, const Var &X2) throw (Error)  {
  return mixture(std::vector<double> {wX1, wX2},
                 std::vector<Var> {X1, X2});
}

Var
stochbb::mixture(double wX1, const Var &X1, double wX2, const Var &X2, double wX3, const Var &X3) throw (Error) {
  return mixture(std::vector<double> {wX1, wX2, wX3},
                 std::vector<Var> {X1, X2, X3});
}

Var
stochbb::mixture(const std::vector<double> &weights, const std::vector<Var> &variables) throw (Error) {
  return new MixtureObj(weights, variables);
}


/* ********************************************************************************************* *
 * Implementation of conditional
 * ********************************************************************************************* */
Var
stochbb::conditional(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2) throw (Error) {
  std::vector<Var> args = {X1, X2};
  splitCommon(args);
  if (! independent(args)) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? Y1 : Y2."
           " Variables X1, X2 are not mutually independent.";
    throw err;
  }
  return new ConditionalObj(args[0], args[1], Y1, Y2);
}


/* ********************************************************************************************* *
 * Implementation of condchain
 * ********************************************************************************************* */
Var
stochbb::condchain(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2) throw (Error) {
  std::vector<Var> args = {X1, X2};
  Var common = splitCommon(args);
  // Check for independence (Y1 and Y2 are allowed to be dependent RVs).
  if (! independent(std::vector<Var> {args[0], args[1], Y1})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? X1+Y1 : X2+Y2."
           " Variables X1, X2, Y1 are not mutually independent.";
    throw err;
  }
  if (! independent(std::vector<Var> {args[0], args[1], Y2})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? X1+Y1 : X2+Y2."
           " Variables X1, X2, Y2 are not mutually independent.";
    throw err;
  }
  return Var(new CondChainObj(args[0], args[1], Y1, Y2)) + common;
}


/* ********************************************************************************************* *
 * Implementation of directConvolve
 * ********************************************************************************************* */
Density
stochbb::directConvolve(const std::vector<Density> &densities) {
  return new ConvolutionDensityObj(densities);
}

Density
stochbb::directConvolve(const Density &a, const Density &b) {
  return directConvolve(std::vector<Density> {a,b});
}

Density
stochbb::directConvolve(const Density &a, const Density &b, const Density &c) {
  return directConvolve(std::vector<Density> {a,b,c});
}

std::ostream &
operator<<(std::ostream &stream, const stochbb::Container &x) {
  if (x.isNull()) {
    stream << "<null>";
  } else {
    x->print(stream);
  }
  return stream;
}

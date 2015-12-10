#include "operators.hh"
#include "randomvariable.hh"
#include "affinetrafo.hh"
#include "chain.hh"
#include "minmax.hh"
#include "mixture.hh"

using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of sbb::delta()
 * ********************************************************************************************* */
Var
stochbb::delta(double value) {
  return AtomicVar::delta(value);
}

/* ********************************************************************************************* *
 * Implementation of sbb::uniform()
 * ********************************************************************************************* */
Var
stochbb::uniform(double a, double b) {
  return AtomicVar::unif(a, b);
}

/* ********************************************************************************************* *
 * Implementation of sbb::normal()
 * ********************************************************************************************* */
Var
stochbb::normal(double mu, double sigma, const std::string &name) {
  return AtomicVar::norm(mu, sigma, name);
}

Var
stochbb::normal(const Var &mu, double sigma, const std::string &name) {
  // If mu is delta distributed -> simplify to atomic random variable
  if (DeltaDensityObj *delta = dynamic_cast<DeltaDensityObj *>(*mu.density())) {
    return stochbb::normal(delta->delay(), sigma, name);
  }
  // Otherwise assemble CompoundVar
  return Compound::norm(mu, AtomicVar::delta(sigma), name);
}

Var
stochbb::normal(double mu, const Var &sigma, const std::string &name) {
  // If sigma is delta distributed -> simplify to atomic random variable
  if (DeltaDensityObj *delta = dynamic_cast<DeltaDensityObj *>(*sigma.density())) {
    return stochbb::normal(mu, delta->delay(), name);
  }
  // Otherwise assemble CompoundVar
  return Compound::norm(AtomicVar::delta(mu), sigma, name);
}

Var
stochbb::normal(const Var &mu, const Var &sigma, const std::string &name) {
  // If mu is delta distributed -> simplify to atomic random variable
  if (DeltaDensityObj *delta = dynamic_cast<DeltaDensityObj *>(*mu.density())) {
    return stochbb::normal(delta->delay(), sigma, name);
  }
  // else if sigma is delta distributed -> simplify to atomic random variable
  if (DeltaDensityObj *delta = dynamic_cast<DeltaDensityObj *>(*sigma.density())) {
    return stochbb::normal(mu, delta->delay(), name);
  }
  // Otherwise assemble CompoundVar
  return Compound::norm(mu, sigma, name);
}


/* ********************************************************************************************* *
 * Implementation of sbb::gamma()
 * ********************************************************************************************* */
Var
stochbb::gamma(double k, double theta, const std::string &name) {
  return AtomicVar::gamma(k, theta, name);
}

Var
stochbb::gamma(const Var &k, double theta, const std::string &name) {
  // If mu is delta distributed -> simplify to atomic random variable
  if (DeltaDensityObj *delta = dynamic_cast<DeltaDensityObj *>(*k.density())) {
    return stochbb::gamma(delta->delay(), theta, name);
  }
  // Otherwise assemble CompoundVar
  return Compound::gamma(k, AtomicVar::delta(theta), name);
}

Var
stochbb::gamma(double k, const Var &theta, const std::string &name) {
  // If sigma is delta distributed -> simplify to atomic random variable
  if (DeltaDensityObj *delta = dynamic_cast<DeltaDensityObj *>(*theta.density())) {
    return stochbb::gamma(k, delta->delay(), name);
  }
  // Otherwise assemble CompoundVar
  return Compound::gamma(AtomicVar::delta(k), theta, name);
}

Var
stochbb::gamma(const Var &k, const Var &theta, const std::string &name) {
  // If mu is delta distributed -> simplify to atomic random variable
  if (DeltaDensityObj *delta = dynamic_cast<DeltaDensityObj *>(*k.density())) {
    return stochbb::gamma(delta->delay(), theta, name);
  }
  // else if sigma is delta distributed -> simplify to atomic random variable
  if (DeltaDensityObj *delta = dynamic_cast<DeltaDensityObj *>(*theta.density())) {
    return stochbb::gamma(k, delta->delay(), name);
  }
  // Otherwise assemble CompoundVar
  return Compound::gamma(k, theta, name);
}


/* ********************************************************************************************* *
 * Implementation of minimum
 * ********************************************************************************************* */
Var stochbb::minimum(const Var &X) {
  std::vector<Var> vars; vars.reserve(1);
  vars.push_back(X);
  return minimum(vars);
}

Var stochbb::minimum(const Var &X1, const Var &X2) {
  std::vector<Var> vars; vars.reserve(2);
  vars.push_back(X1); vars.push_back(X2);
  return minimum(vars);
}

Var stochbb::minimum(const Var &X1, const Var &X2, const Var &X3) {
  std::vector<Var> vars; vars.reserve(3);
  vars.push_back(X1); vars.push_back(X2); vars.push_back(X3);
  return minimum(vars);
}

Var
stochbb::minimum(const std::vector<Var> &variables) {
  // Check size of variables
  if (0 == variables.size()) {
    AssumptionError err;
    err << "Cannot construct the minimum of no variable.";
    throw err;
  } else if (1 == variables.size()) {
    // If only one variable is given -> return it
    return variables[0];
  }

  // Expand minimum objects, e.g. min(min(A,B),C) -> min(A,B,C)
  std::vector<Var> vars; vars.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    if (variables[i].is<Minimum>()) {
      Minimum max = variables[i].as<Minimum>();
      for (size_t j=0; j<max.numVariables(); j++) {
        vars.push_back(max.variable(j));
      }
    }
  }

  // Find the common part of the variables formed as a sum
  std::vector<VarSet> indepvars;
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
  VarObj *result = 0;
  std::vector<Var> args; args.reserve(indepvars.size());
  for (size_t i=0; i<indepvars.size(); i++) {
    if (1 == indepvars[i].size()) {
      args.push_back(*indepvars[i].begin());
    } else {
      std::vector<Var> chain_args; chain_args.reserve(indepvars[i].size());
      VarSet::iterator item = indepvars[i].begin();
      for (; item != indepvars[i].end(); item++) {
        chain_args.push_back(*item);
      }
      args.push_back(new ChainObj(chain_args));
    }
    result = new MinimumObj(args);
  }

  // If there is a common part -> assemble chain
  if (! common.isEmpty()) {
    std::vector<Var> chain_args; chain_args.reserve(common.size()+1);
    chain_args.push_back(result);
    VarSet::iterator item = common.begin();
    for (; item != common.end(); item++) {
      chain_args.push_back(*item);
    }
    result = new ChainObj(chain_args);
  }

  // done.
  return result;
}


/* ********************************************************************************************* *
 * Implementation of maximum
 * ********************************************************************************************* */
Var stochbb::maximum(const Var &X) {
  std::vector<Var> vars; vars.reserve(1);
  vars.push_back(X);
  return maximum(vars);
}

Var stochbb::maximum(const Var &X1, const Var &X2) {
  std::vector<Var> vars; vars.reserve(2);
  vars.push_back(X1); vars.push_back(X2);
  return maximum(vars);
}

Var stochbb::maximum(const Var &X1, const Var &X2, const Var &X3) {
  std::vector<Var> vars; vars.reserve(3);
  vars.push_back(X1); vars.push_back(X2); vars.push_back(X3);
  return maximum(vars);
}

Var
stochbb::maximum(const std::vector<Var> &variables) {
  // Check size of variables
  if (0 == variables.size()) {
    AssumptionError err;
    err << "Cannot construct the maximum of no variable.";
    throw err;
  }

  // If only one variable is given -> return it
  if (1 == variables.size()) {
    return variables[0];
  }

  // Expand maximum objects, e.g. max(max(A,B),C) -> max(A,B,C)
  std::vector<Var> vars; vars.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    if (variables[i].is<Maximum>()) {
      Maximum max = variables[i].as<Maximum>();
      for (size_t j=0; j<max.numVariables(); j++) {
        vars.push_back(max.variable(j));
      }
    }
  }

  // Find the common part of the variables formed as a sum
  std::vector<VarSet> indepvars;
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
  VarObj *result = 0;
  std::vector<Var> args; args.reserve(indepvars.size());
  for (size_t i=0; i<indepvars.size(); i++) {
    if (1 == indepvars[i].size()) {
      args.push_back(*indepvars[i].begin());
    } else {
      std::vector<Var> chain_args; chain_args.reserve(indepvars[i].size());
      VarSet::iterator item = indepvars[i].begin();
      for (; item != indepvars[i].end(); item++) {
        chain_args.push_back(*item);
      }
      args.push_back(new ChainObj(chain_args));
    }
    result = new MaximumObj(args);
  }

  // If there is a common part -> assemble chain
  if (! common.isEmpty()) {
    std::vector<Var> chain_args; chain_args.reserve(common.size()+1);
    chain_args.push_back(result);
    VarSet::iterator item = common.begin();
    for (; item != common.end(); item++) {
      chain_args.push_back(*item);
    }
    result = new ChainObj(chain_args);
  }

  // done.
  return result;
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


/* ********************************************************************************************* *
 * Implementation of chain
 * ********************************************************************************************* */
Var
stochbb::chain(const std::vector<Var> &vars) {
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
stochbb::mixture(double wX1, const Var &X1, double wX2, const Var &X2) {
  std::vector<double> weights; weights.reserve(2);
  std::vector<Var> variables; variables.reserve(2);
  weights.push_back(wX1); weights.push_back(wX2);
  variables.push_back(X1); variables.push_back(X2);
  return mixture(weights, variables);
}

Var
stochbb::mixture(double wX1, const Var &X1, double wX2, const Var &X2, double wX3, const Var &X3) {
  std::vector<double> weights; weights.reserve(3);
  std::vector<Var> variables; variables.reserve(3);
  weights.push_back(wX1); weights.push_back(wX2); weights.push_back(wX3);
  variables.push_back(X1); variables.push_back(X2); variables.push_back(X3);
  return mixture(weights, variables);
}

Var
stochbb::mixture(const std::vector<double> &weights, const std::vector<Var> &variables) {
  return new MixtureObj(weights, variables);
}

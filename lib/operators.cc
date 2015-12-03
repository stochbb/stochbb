#include "operators.hh"
#include "minmax.hh"
#include "chain.hh"

using namespace sbb;


/* ********************************************************************************************* *
 * Implementation of minimum
 * ********************************************************************************************* */
Var
sbb::minimum(const std::vector<Var> &variables) {
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
Var
sbb::maximum(const std::vector<Var> &variables) {
  // Check size of variables
  if (0 == variables.size()) {
    AssumptionError err;
    err << "Cannot construct the maximum of no variable.";
    throw err;
  } else if (1 == variables.size()) {
    // If only one variable is given -> return it
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



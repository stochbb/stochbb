#ifndef __SBB_COMPOUND_HH__
#define __SBB_COMPOUND_HH__

#include "api.hh"
#include "randomvariable.hh"
#include "distribution.hh"
#include <Eigen/Eigen>


namespace stochbb {

/** Base class of all compound random variables. */
class CompoundObj : public DerivedVarObj
{
protected:
  /** Hidden constructor.
   * @param vars Specifies the variables, the compound depends on.
   * @param name Specifies the optional name for the random variable. */
  CompoundObj(const std::vector<Var> &vars, const std::string &name);

public:
  virtual void mark();
};

/** Implements the generic density of compound random variables. That is a random variable
 * that is distributed with an atomic density where the parameters are random variables too. */
class GenericCompoundDensityObj: public DensityObj
{
public:
  /** Constructs a compound density with the given distribution family and parameter distributions. */
  GenericCompoundDensityObj(DistributionObj *dist, const std::vector<DensityObj *> &params);
  /** Destructor. */
  virtual ~GenericCompoundDensityObj();
  void mark();

  /** Returns a reference to the distribution object. */
  Distribution distribution() const;
  /** Returns the number of parameters of the density. */
  size_t nParams() const;
  /** Returns a specific parameter density. */
  Density parameter(size_t i) const;

  void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;

  Density affine(double scale, double shift) const;

  void rangeEst(double alpha, double &a, double &b) const;

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

protected:
  /** Helper function to iterate over several indices.
   * @param i Specifies the running index.
   * @param N Specifies the number of dimensions per index.
   * @param idxs On exit, contains the vector of indices. */
  void _to_param_indices(size_t i, size_t N, std::vector<size_t> &idxs) const;

protected:
  /** Holds a reference to the distribution instance. */
  DistributionObj *_distribution;
  /** Holds references to the parameter densities. */
  std::vector<DensityObj *> _parameters;
};


/** Base class of all compound random variables. */
class GenericCompoundObj : public CompoundObj
{
public:
  /** Constructs a new compound random variable object.
   * @param vars Specifies the variables, the compound depends on.
   * @param name Specifies the optional name for the random variable. */
  GenericCompoundObj(const std::vector<Var> &vars, const Distribution &distribution, const std::string &name);

  virtual void mark();

  virtual Density density();
  Distribution distribution();

protected:
  /** A reference to the density. */
  GenericCompoundDensityObj *_density;
};

}

#endif // COMPOUND_HH

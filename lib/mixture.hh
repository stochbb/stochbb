#ifndef __SBB_MIXTURE_HH__
#define __SBB_MIXTURE_HH__

#include "api.hh"
#include "randomvariable.hh"
#include "density.hh"


namespace stochbb {

/** Implements a mixture density.
 * @ingroup density */
class MixtureDensityObj: public DensityObj
{
protected:
  /** Hidden constructor from the given weights and densities. */
  MixtureDensityObj(const std::vector<double> &weights, const std::vector<DensityObj *> &densities,
                    double scale=1, double shift=0);

public:
  /** Constructs a mixture density from the given weights and variables. */
  MixtureDensityObj(const std::vector<double> &weights, const std::vector<VarObj *> &variables,
                    double scale=1, double shift=0);

  virtual void mark();
  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;
  virtual void rangeEst(double alpha, double &a, double &b) const;
  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

protected:
  /** The vector of weights (normalized). */
  std::vector<double> _weights;
  /** The vector of densities of the random variables being mixed. */
  std::vector<DensityObj *> _densities;
  /** The scale of the affine transform. */
  double _scale;
  /** The shift of the affine transform. */
  double _shift;
};


/** Represents a random variable being a mixture of other random variables.
 * @ingroup rv */
class MixtureObj : public DerivedVarObj
{
public:
  /** Constructs a mixture of the given @c variables and @c weights. */
  MixtureObj(const std::vector<double> &weights, const std::vector<Var> &variables,
             const std::string &name="") throw (Error);

  virtual void mark();
  virtual Density density();

  /** Returns the weight of the i-th random variable. */
  inline double weight(size_t i) const { return _weights[i]; }

  virtual void sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                      Eigen::Ref<Eigen::MatrixXd> samples) const;
  virtual void print(std::ostream &stream) const;

protected:
  /** The vector of weights. */
  std::vector<double> _weights;
  /** The density of the mixture. */
  MixtureDensityObj *_density;
};

}

#endif // __SBB_MIXTURE_HH__

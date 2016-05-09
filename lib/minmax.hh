#ifndef __SBB_MINMAX_HH__
#define __SBB_MINMAX_HH__

#include <vector>
#include "api.hh"
#include "randomvariable.hh"

namespace stochbb {

/** Implements the density of a random variable being the maximum of several independent random
 * variables.
 * @ingroup density */
class MaximumDensityObj: public DensityObj
{
protected:
  /** Hidden constructor.
   * @param densities Specifies the vector of densities associated with the independent random
   *   variables (weak references).
   * @param scale Specifies the optional scaling of the density. Default @c scale = 1.
   * @param shift Specifies the optional shift of the density. Default @c shift = 0.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MaximumDensityObj(const std::vector<DensityObj *> &densities, double scale=1, double shift=0);

  /** Hidden constructor.
   * @param variables Specifies the vector of independent random variables (weak references).
   * @param scale Specifies the optional scaling of the density. Default @c scale = 1.
   * @param shift Specifies the optional shift of the density. Default @c shift = 0.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MaximumDensityObj(const std::vector<VarObj *> &variables, double scale=1, double shift=0);

public:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables.
   * @param scale Specifies the optional scaling of the density. Default @c scale = 1.
   * @param shift Specifies the optional shift of the density. Default @c shift = 0.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MaximumDensityObj(const std::vector<Var> &variables, double scale=1, double shift=0);

  /** Destructor. */
  virtual ~MaximumDensityObj();

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;
  virtual void rangeEst(double alpha, double &a, double &b) const;

  /** Returns a vector of weak references to the underlaying densities. */
  inline const std::vector<DensityObj *> &densities() const { return _densities; }

  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

protected:
  /** The vector of densities. */
  std::vector<DensityObj *> _densities;
  /** Scale of the affine transform. */
  double _scale;
  /** Shift of the affine transform. */
  double _shift;

  friend class MaximumObj;
};


/** Implements the density of a random variable being the minimum of several independent random
 * variables.
 * @ingroup density */
class MinimumDensityObj: public DensityObj
{
protected:
  /** Hidden constructor.
   * @param densities Specifies the vector of densities associated with the independent random
   *   variables (weak references).
   * @param scale Specifies the optional scaling of the density. Default @c scale = 1.
   * @param shift Specifies the optional shift of the density. Default @c shift = 0.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MinimumDensityObj(const std::vector<DensityObj *> &densities, double scale=1, double shift=0);
  /** Hidden constructor.
   * @param variables Specifies the vector of independent random variables (weak references).
   * @param scale Specifies the optional scaling of the density. Default @c scale = 1.
   * @param shift Specifies the optional shift of the density. Default @c shift = 0.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MinimumDensityObj(const std::vector<VarObj *> &variables, double scale=1, double shift=0);

public:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables.
   * @param scale Specifies the optional scaling of the density. Default @c scale = 1.
   * @param shift Specifies the optional shift of the density. Default @c shift = 0
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MinimumDensityObj(const std::vector<Var> &variables, double scale=1, double shift=0);

  /** Destructor. */
  virtual ~MinimumDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;
  virtual void rangeEst(double alpha, double &a, double &b) const;

  /** Returns a vector of weak references to the underlaying densities. */
  inline const std::vector<DensityObj *> &densities() const { return _densities; }

  virtual int compare(const DensityObj &other) const;
  virtual void print(std::ostream &stream) const;

protected:
  /** The vector of densities. */
  std::vector<DensityObj *> _densities;
  /** Scale of the affine transform. */
  double _scale;
  /** Shift of the affine transform. */
  double _shift;
  friend class MinimumObj;
};


/** Implements the random variable being the maximum of several independent random
 * variables.
 * @ingroup rv */
class MaximumObj: public DerivedVarObj
{
public:
  /** Constructor from a vector of random variables.
   * @throws AssumptionError If these random variables are not independent. */
  MaximumObj(const std::vector<Var> &variables, const std::string &name="");

  /** Destructor. */
  virtual ~MaximumObj();

  virtual void mark();

  virtual Density density();

  virtual void sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                      Eigen::Ref<Eigen::MatrixXd> samples) const;

  virtual void print(std::ostream &stream) const;

protected:
  /** The density object. */
  MaximumDensityObj *_density;
};


/** Implements the random variable being the minimum of several independent random
 * variables.
 * @ingroup rv */
class MinimumObj: public DerivedVarObj
{
public:
  /** Constructor from a vector of random variables.
   * @throws AssumptionError If these random variables are not independent. */
  MinimumObj(const std::vector<Var> &variables, const std::string &name="");

  /** Destructor. */
  virtual ~MinimumObj();

  virtual void mark();

  virtual Density density();

  virtual void sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                      Eigen::Ref<Eigen::MatrixXd> samples) const;

  virtual void print(std::ostream &stream) const;

protected:
  /** The density object of this random variable. */
  MinimumDensityObj *_density;
};

}

#endif // __SBB_MINMAX_HH__

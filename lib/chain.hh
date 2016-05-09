#ifndef CHAIN_HH
#define CHAIN_HH

#include <vector>
#include "api.hh"
#include "randomvariable.hh"

namespace stochbb {

/** Implements the convolution of several PDFs using the FFT convolution.
 * @ingroup density */
class ConvolutionDensityObj: public DensityObj
{
protected:
  /** Constructs a new PDF as the convolution of the given PDFs (weak references). */
  ConvolutionDensityObj(const std::vector<DensityObj *> &densities, double scale=1, double shift=0);

public:
  /** Constructs a new PDF as the convolution of the PDFs of the given variables. */
  ConvolutionDensityObj(const std::vector<Var> &variables, double scale=1, double shift=0);
  /** Constructs a new PDF as the convolution of the given PDFs. */
  ConvolutionDensityObj(const std::vector<Density> &densities, double scale=1, double shift=0);
  /** Destructor. */
  virtual ~ConvolutionDensityObj();

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const;
  virtual Density affine(double scale, double shift) const;
  virtual void rangeEst(double alpha, double &a, double &b) const;

  /** Returns the number of underlaying densities. */
  inline size_t numDensities() const { return _densities.size(); }

  /** Returns the i-th density. */
  inline Density density(size_t i) const {
    _densities[i]->ref();
    return _densities[i];
  }

  /** Returns a vector of weak references to the densities of the underlaying variables. */
  inline const std::vector<DensityObj *> &densities() const {
    return _densities;
  }

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

protected:
  /** The @c DensityObj instances of the underlaying variables. */
  std::vector<DensityObj *> _densities;
  /** Scale of the affine transform. */
  double _scale;
  /** Shift of the affine transform. */
  double _shift;

  friend class ChainObj;
};


/** Returns the convoultion of the given densities.
 * This function first tries to perform the convolutions analytically and
 * resorts to the numerical convolution if necessary.
 * @ingroup density */
Density convolve(const std::vector<Density> &densities, double scale=1, double shift=0);

/** Represetns the sum of several independent random variables.
 * @ingroup rv */
class ChainObj : public DerivedVarObj
{
public:
  /** Constructs the sum of the given random variables. */
  ChainObj(const std::vector<Var> &variables, const std::string &name="");
  /** Destructor. */
  virtual ~ChainObj();
  virtual void mark();

  virtual Density density();

  virtual void print(std::ostream &stream) const;
  virtual void sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                      Eigen::Ref<Eigen::MatrixXd> samples) const;

protected:
  /** The density of the sum, the convolution of all PDFs of the underlaying random variables. */
  DensityObj *_density;
};

}

#endif // CHAIN_HH

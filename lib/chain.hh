#ifndef CHAIN_HH
#define CHAIN_HH

#include <vector>
#include "api.hh"
#include "randomvariable.hh"

namespace sbb {

/** Implements the convolution of several PDFs using the FFT convolution.
 * @ingroup internal */
class ConvolutionDensityObj: public DensityObj
{
public:
  /** Constructs a new PDF as the convolution of the PDFs of the given variables. */
  ConvolutionDensityObj(const std::vector<VarObj *> &variables);
  /** Destructor. */
  virtual ~ConvolutionDensityObj();

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

  /** Returns the number of underlaying densities. */
  inline size_t numDensities() const { return _densities.size(); }
  /** Returns the i-th density. */
  inline Density density(size_t i) const {
    _densities[i]->ref();
    return _densities[i];
  }

  /** Returns the densities of the underlaying variables. */
  inline const std::vector<DensityObj *> &densities() const {
    return _densities;
  }

protected:
  /** The @c DensityObj instances of the underlaying variables. */
  std::vector<DensityObj *> _densities;
};


/** Represetns the sum of several independent random variables.
 * @ingroup internal */
class ChainObj : public VarObj
{
public:
  /** Constructs the sum of the given random variables. */
  ChainObj(VarObj *a, VarObj *b, const std::string &name="");
  /** Constructs the sum of the given random variables. */
  ChainObj(const std::vector<VarObj *> &variables, const std::string &name="");
  /** Destructor. */
  virtual ~ChainObj();
  virtual void mark();

  virtual DensityObj *density();

  inline size_t numVariables() const { return _variables.size(); }
  inline Var variable(size_t i) { return _variables[i]; }

protected:
  /** References to the underlaying random variables. */
  std::vector<VarObj *> _variables;
  /** The density of the sum, the convolution of all PDFs of the underlaying random variables. */
  ConvolutionDensityObj *_density;
};

}

#endif // CHAIN_HH

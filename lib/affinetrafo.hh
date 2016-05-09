#ifndef __SBB_AFFINETRAFO_HH__
#define __SBB_AFFINETRAFO_HH__

#include "api.hh"
#include "randomvariable.hh"

namespace stochbb {

/** Represents an affine transformation of a random variable.
 * That is, \f$Y = a\,X + b\f$, where the real parameter \f$a\f$ is called the scale (dillatation)
 * and \f$b\f$ the shift (translation). If \f$X\sim f(x)\f$ then
 * \f$Y = f\left(\frac{x-b}{a}\right)\f$.
 * @ingroup rv */
class AffineTrafoObj : public DerivedVarObj
{
public:
  /** Constructs a affine transformation of the given @c variable with the given @c scale and
   * @c shift. */
  AffineTrafoObj(double scale, double shift, const Var &variable, const std::string &name="");

  virtual void mark();

  virtual Density density();

  /** Returns the scale factor of the affine transform. */
  inline double scale() const { return _scale; }
  /** Returns the shift of the affine transform. */
  inline double shift() const { return _shift; }

  virtual void print(std::ostream &stream) const;

  virtual void sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                      Eigen::Ref<Eigen::MatrixXd> samples) const;

protected:
  /** The density of the scaled and shifted RV. */
  DensityObj *_density;
  /** The scale parameter. */
  double _scale;
  /** The shift parameter. */
  double _shift;
};

}

#endif // AFFINETRAFO_HH

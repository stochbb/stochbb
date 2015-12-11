#ifndef __SBB_AFFINETRAFO_HH__
#define __SBB_AFFINETRAFO_HH__

#include "api.hh"
#include "randomvariable.hh"

namespace stochbb {

class AffineTrafoObj : public DerivedVarObj
{
public:
  AffineTrafoObj(double scale, double shift, const Var &variable, const std::string &name="");

  virtual void mark();

  virtual Density density();

  /** Returns the scale factor of the affine transform. */
  inline double scale() const { return _scale; }
  /** Returns the shift of the affine transform. */
  inline double shift() const { return _shift; }

  virtual void print(std::ostream &stream) const;

protected:
  DensityObj *_density;
  double _scale;
  double _shift;
};

}

#endif // AFFINETRAFO_HH

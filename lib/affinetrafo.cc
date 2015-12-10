#include "affinetrafo.hh"

using namespace stochbb;


AffineTrafoObj::AffineTrafoObj(double scale, double shift, const Var &variable, const std::string &name)
  : DerivedVarObj(std::vector<Var>(1,variable), name), _density(0), _scale(scale), _shift(shift)
{
  _density = *variable->density().affine(scale,shift);
}

void
AffineTrafoObj::mark() {
  if (isMarked()) { return; }
  DerivedVarObj::mark();
  if (_density) { _density->mark(); }
}

Density
AffineTrafoObj::density() {
  _density->ref();
  return _density;
}

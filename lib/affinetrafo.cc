#include "affinetrafo.hh"

using namespace stochbb;


AffineTrafoObj::AffineTrafoObj(double scale, double shift, const Var &variable, const std::string &name)
  : DerivedVarObj(std::vector<Var> { variable }, name), _density(0), _scale(scale), _shift(shift)
{
  _density = *variable->density().affine(scale, shift);
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

void
AffineTrafoObj::print(std::ostream &stream) const {
  stream << "<AffineTrafo " << _scale << " * ";
  _variables[0]->print(stream);
  stream << " + " << _shift << " #" << this << ">";
}

void
AffineTrafoObj::sample(size_t outIdx, const Eigen::Ref<IndexVector> &indices,
                       Eigen::Ref<Eigen::MatrixXd> samples) const
{
  samples.col(outIdx) = _scale*samples.col(indices(0)).array()+_shift;
}

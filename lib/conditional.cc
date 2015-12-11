#include "conditional.hh"
#include "operators.hh"

using namespace stochbb;


/* ********************************************************************************************* *
 * Implementation of ConditionalDensityObj
 * ********************************************************************************************* */
ConditionalDensityObj::ConditionalDensityObj(DensityObj *X1, DensityObj *X2, DensityObj *Y1, DensityObj *Y2)
  : DensityObj(), _X1(X1), _X2(X2), _Y1(Y1), _Y2(Y2)
{
  // pass...
}

ConditionalDensityObj::ConditionalDensityObj(const Var& X1, const Var &X2, const Var &Y1, const Var &Y2)
  : DensityObj(), _X1(*X1.density()), _X2(*X2.density()), _Y1(*Y1.density()), _Y2(*Y2.density())
{
  // pass...
}

void
ConditionalDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  _X1->mark(); _X2->mark();
  _Y1->mark(); _Y2->mark();
}

void
ConditionalDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  Eigen::VectorXd f1(out.size()), f2(out.size());
  _X1->evalCDF(Tmin, Tmax, f1);
  _X2->eval(Tmin, Tmax, f2);
  double dt = (Tmax-Tmin)/out.size();
  double p = dt*((f1.array()*f2.array()).sum());
  _Y1->eval(Tmin, Tmax, f1);
  _Y2->eval(Tmin, Tmax, f2);
  out = p*f1 + (1-p)*f2;
}

void
ConditionalDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  Eigen::VectorXd f1(out.size()), f2(out.size());
  _X1->evalCDF(Tmin, Tmax, f1);
  _X2->eval(Tmin, Tmax, f2);
  double dt = (Tmax-Tmin)/out.size();
  double p = dt*((f1.array()*f2.array()).sum());
  _Y1->evalCDF(Tmin, Tmax, f1);
  _Y2->evalCDF(Tmin, Tmax, f2);
  out = p*f1 + (1-p)*f2;
}

Density
ConditionalDensityObj::affine(double scale, double shift) const {
  // Apply affine transformation only on the Y* densities, the condition
  // RV densities are unaffected
  return new ConditionalDensityObj(
        _X1, _X2, *_Y1->affine(scale, shift), *_Y2->affine(scale, shift));
}


/* ********************************************************************************************* *
 * Implementation of ConditionalObj
 * ********************************************************************************************* */
ConditionalObj::ConditionalObj(const Var& X1, const Var &X2, const Var &Y1, const Var &Y2,
                               const std::string &name)
  : DerivedVarObj(std::vector<Var> {X1, X2, Y1, Y2}, name), _density(0)
{
  // Check for independence (Y1 and Y2 are allowed to be dependent RVs).
  //if (! independent(std::vector<Var> {X1, X2, Y1})) {
  if (! independent(std::vector<Var> {X1, X2})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? Y1 : Y2."
           " Variables X1, X2, Y1 are not mutually independent.";
    throw err;
  }
  /*if (! independent(std::vector<Var> {X1, X2, Y2})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? Y1 : Y2."
           " Variables X1, X2, Y2 are not mutually independent.";
    throw err;
  }*/

  _density = new ConditionalDensityObj(X1, X2, Y1, Y2);
}


void
ConditionalObj::mark() {
  if (isMarked()) { return; }
  DerivedVarObj::mark();
  if (_density) { _density->mark(); }
}

Density
ConditionalObj::density() {
  _density->ref();
  return _density;
}



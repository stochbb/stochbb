#include "conditional.hh"
#include "operators.hh"
#include <unsupported/Eigen/FFT>

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
  : DensityObj(), _X1(0), _X2(0), _Y1(0), _Y2(0)
{
  _X1 = *X1.density(); _X2 = *X2.density();
  _Y1 = *Y1.density(); _Y2 = *Y2.density();
}

void
ConditionalDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  if (_X1) { _X1->mark(); }
  if (_X2) { _X2->mark(); }
  if (_Y1) { _Y1->mark(); }
  if (_Y2) { _Y2->mark(); }
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
  if (! independent(std::vector<Var> {X1, X2, Y1})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? Y1 : Y2."
           " Variables X1, X2, Y1 are not mutually independent.";
    throw err;
  }
  if (! independent(std::vector<Var> {X1, X2, Y2})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? Y1 : Y2."
           " Variables X1, X2, Y2 are not mutually independent.";
    throw err;
  }

  _density = new ConditionalDensityObj(X1, X2, Y1, Y2);
  _density->unref();
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


/* ********************************************************************************************* *
 * Implementation of CondChainDensityObj
 * ********************************************************************************************* */
CondChainDensityObj::CondChainDensityObj(DensityObj *X1, DensityObj *X2, DensityObj *Y1, DensityObj *Y2)
  : DensityObj(), _X1(X1), _X2(X2), _Y1(Y1), _Y2(Y2)
{
  // pass...
}

CondChainDensityObj::CondChainDensityObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2)
  : DensityObj(), _X1(0), _X2(0), _Y1(0), _Y2(0)
{
  _X1 = *X1.density(); _X2 = *X2.density();
  _Y1 = *Y1.density(); _Y2 = *Y2.density();
}

void
CondChainDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  if (_X1) { _X1->mark(); }
  if (_X2) { _X2->mark(); }
  if (_Y1) { _Y1->mark(); }
  if (_Y2) { _Y2->mark(); }
}

void
CondChainDensityObj::eval(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  Eigen::FFT<double> fft;
  Eigen::VectorXd tmp1(2*out.size());  tmp1.setZero();
  Eigen::VectorXcd tmp2(2*out.size());
  Eigen::VectorXcd prod(2*out.size()); //prod.setOnes();
  Eigen::VectorXd sum(out.size());

  /// @bug Fix shift if Tmin != 0, like in @c ConvolutionDensity!
  logDebug() << "Eval CondChainDensityObj on int. [" << Tmin << ", " << Tmax << "].";

  double dt = (Tmax-Tmin)/out.size();
  // Get sample rate
  double Fs = out.size()/(Tmax-Tmin);
  // Compute time-shift
  std::complex<double> dphi(0, 2*M_PI*Tmin*Fs/double(2*out.size()));

  // Perform FFT convolution
  // compute FFT(f_{X_1}*F_{X_2})
  _X1->eval(Tmin, Tmax, out);
  tmp1.head(out.size()) = out*dt;
  _X2->evalCDF(Tmin, Tmax, out);
  tmp1.head(out.size()).array() *= (1-out.array());
  fft.fwd(prod, tmp1);
  // apply time shift
  for (int i=1; i<out.size();i++) {
    prod[i] *= std::exp(-double(i)*dphi);
    prod[2*out.size()-i] *= std::exp(double(i)*dphi);
  }
  // compute FFT(f_{Y_1})
  _Y1->eval(Tmin, Tmax, out);
  tmp1.head(out.size()) = out*dt;
  fft.fwd(tmp2, tmp1);
  // apply time shift
  for (int i=1; i<out.size();i++) {
    tmp2[i] *= std::exp(-double(i)*dphi);
    tmp2[2*out.size()-i] *= std::exp(double(i)*dphi);
  }
  // FFT_INV( FFT(f_{X_1}*F_{X_2}) * FFT(f_{Y_1}) )
  prod.array() *= tmp2.array();
  // Reverse time-shift
  for (int i=1; i<out.size();i++) {
    prod[i] *= std::exp(double(i)*dphi);
    prod[2*out.size()-i] *= std::exp(-double(i)*dphi);
  }
  fft.inv(tmp1, prod);
  sum.noalias() = tmp1.head(out.size())/dt;

  tmp1.setZero();
  // compute FFT(f_{X_2}*F_{X_1})
  _X2->eval(Tmin, Tmax, out);
  tmp1.head(out.size()) = out*dt;
  _X1->evalCDF(Tmin, Tmax, out);
  tmp1.head(out.size()).array() *= (1-out.array());
  fft.fwd(prod, tmp1);
  // apply time shift
  for (int i=1; i<out.size();i++) {
    prod[i] *= std::exp(-double(i)*dphi);
    prod[2*out.size()-i] *= std::exp(double(i)*dphi);
  }
  // compute FFT(f_{Y_2})
  _Y2->eval(Tmin, Tmax, out);
  tmp1.head(out.size()) = out*dt;
  fft.fwd(tmp2, tmp1);
  // apply time shift
  for (int i=1; i<out.size();i++) {
    tmp2[i] *= std::exp(-double(i)*dphi);
    tmp2[2*out.size()-i] *= std::exp(double(i)*dphi);
  }
  // FFT_INV( FFT(f_{X_2}*F_{X_1}) * FFT(f_{Y_2}) )
  prod.array() *= tmp2.array();
  // Reverse time-shift
  for (int i=1; i<out.size();i++) {
    prod[i] *= std::exp(double(i)*dphi);
    prod[2*out.size()-i] *= std::exp(-double(i)*dphi);
  }
  fft.inv(tmp1, prod);
  // sum
  out.noalias() = sum + tmp1.head(out.size())/dt;
}

void
CondChainDensityObj::evalCDF(double Tmin, double Tmax, Eigen::Ref<Eigen::VectorXd> out) const {
  if (0 == out.size()) { return; }
  double dt = (Tmax-Tmin)/(out.size());
  this->eval(Tmin, Tmax, out);
  out[0] *=dt;
  for (int i=1; i<out.size(); i++) {
    out[i] = out[i-1] + out[i]*dt;
  }
}

Density
CondChainDensityObj::affine(double scale, double shift) const {
  // Simply apply affine transform on all densities
  return new CondChainDensityObj(*_X1->affine(scale, 0), *_X2->affine(scale, 0),
                                 *_Y1->affine(scale, shift), *_Y2->affine(scale, shift));
}

void
CondChainDensityObj::print(std::ostream &stream) const {
  stream << "<CondChainDensity of X1=";
  _X1->print(stream);
  stream << ", X2=";
  _X2->print(stream);
  stream << ", Y1=";
  _Y1->print(stream);
  stream << ", Y2=";
  _Y2->print(stream);
  stream << " #" << this << ">";
}


/* ********************************************************************************************* *
 * Implementation of CondChainObj
 * ********************************************************************************************* */
CondChainObj::CondChainObj(const Var &X1, const Var &X2, const Var &Y1, const Var &Y2, const std::string &name)
  : DerivedVarObj(std::vector<Var> {X1, X2, Y1, Y2}, name), _density(0)
{
  // Check for independence (Y1 and Y2 are allowed to be dependent RVs).
  if (! independent(std::vector<Var> {X1, X2, Y1})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? X1+Y1 : X2+Y2."
           " Variables X1, X2, Y1 are not mutually independent.";
    throw err;
  }
  if (! independent(std::vector<Var> {X1, X2, Y2})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? X1+Y1 : X2+Y2."
           " Variables X1, X2, Y2 are not mutually independent.";
    throw err;
  }

  _density = new CondChainDensityObj(X1, X2, Y1, Y2);
  _density->unref();
}

void
CondChainObj::mark() {
  if (isMarked()) { return; }
  DerivedVarObj::mark();
  if (_density) { _density->mark(); }
}

Density
CondChainObj::density() {
  _density->ref();
  return _density;
}

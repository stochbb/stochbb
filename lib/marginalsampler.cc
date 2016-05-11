#include "marginalsampler.hh"
#include "rng.hh"
#include "randomvariable.hh"


using namespace stochbb;

// log(N) binary search algorithm for the index at which p must be inserted to maintain order
inline size_t _find_index(double p, size_t a, size_t b, const Eigen::VectorXd &cdf) {
  while (1 < (b-a)) {
    size_t mid = a+(b-a)/2;
    if (p < cdf[mid]) { b = mid; }
    else { a = mid; }
  }
  return b;
}

MarginalSamplerObj::MarginalSamplerObj(const Var &variable, double Tmin, double Tmax, size_t steps)
  : Object(), _variable(*variable), _Tmin(Tmin), _Tmax(Tmax)
{
  // Resize cdf vector
  _cdf.resize(steps);
  // Evaluate CDF
  _variable->density().evalCDF(_Tmin, _Tmax, _cdf);
}

void
MarginalSamplerObj::mark() {
  if (isMarked()) { return; }
  Object::mark();
  _variable->mark();
}

void
MarginalSamplerObj::sample(Eigen::Ref<Eigen::VectorXd> out) {
  double dt = (_Tmax-_Tmin)/_cdf.size();
  RNG &rng = RNG::get();
  std::uniform_real_distribution<double> sampler(0,1);
  // First, sample probabilities from unif(0,1)
  for (int i=0; i<out.size(); i++) {
    double p = sampler(rng);
    if (p < _cdf[0]) {
      out[i] = _Tmin;
    } else if (p >= _cdf[_cdf.size()-1]) {
      out[i] = _Tmax;
    } else {
      size_t idx = _find_index(p, 0, _cdf.size()-1, _cdf);
      // linear interpolation
      out[i] = _Tmin+dt*(idx + (_cdf[idx]-_cdf[idx-1])/(p-_cdf[idx-1]) - 1);
    }
  }
}


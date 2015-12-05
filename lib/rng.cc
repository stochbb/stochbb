#include "rng.hh"
#include <ctime>
#include <cmath>


using namespace sbb;

RNG *RNG::_instance = 0;

RNG::RNG()
  : MersenneTwister(time(0))
{
  _instance = this;
}

RNG *
RNG::get() {
  if (0 == RNG::_instance) { RNG::_instance = new RNG(); }
  return RNG::_instance;
}

double
RNG::unif() {
  return get()->rand();
}

double
RNG::unif(double a, double b) {
  return a+unif()*(b-a);
}

double
RNG::norm() {
  double u, v, s;
  do {
    u = unif(-1,1);
    v = unif(-1,1);
    s = u*u + v*v;
  } while (s>=1);
  return u*std::sqrt(-2*std::log(s)/s);
}

double
RNG::norm(double mu, double sigma) {
  return mu + sigma*norm();
}

double
RNG::gamma(double a, double b) {
  /* assume a > 0 */
  if (a < 1) {
    double u = unif();
    return gamma(1.0 + a, b) * pow (u, 1.0 / a);
  }

  double x, v, u;
  double d = a - 1.0 / 3.0;
  double c = (1.0 / 3.0) / sqrt (d);

  while (1) {
    do {
      x = norm();
      v = 1.0 + c * x;
    } while (v <= 0);

    v = v * v * v;
    u = unif();

    if (u < (1 - 0.0331 * x * x * x * x)) break;
    if (log (u) < (0.5 * x * x + d * (1 - v + log (v)))) break;
  }

  return b * d * v;
}


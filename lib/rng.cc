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
  double u = unif(-1,1), v = unif(-1,1);
  double s = u*u + v*v;
  return u*std::sqrt(-2*std::log(s)/s);
}

double
RNG::norm(double mean, double stddev) {
  return mean + stddev*norm();
}


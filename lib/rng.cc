#include "rng.hh"
#include <ctime>
#include <cmath>


using namespace stochbb;

RNG *RNG::_instance = 0;

RNG::RNG()
  : std::mt19937_64(time(NULL))
{
  _instance = this;
}

RNG &
RNG::get() {
  if (0 == RNG::_instance) { RNG::_instance = new RNG(); }
  return *RNG::_instance;
}

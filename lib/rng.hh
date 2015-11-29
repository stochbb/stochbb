#ifndef __SBB_RNG_HH__
#define __SBB_RNG_HH__

#include "mersennetwister.hh"

namespace sbb {

class RNG: protected MersenneTwister
{
protected:
  RNG();

public:
  static double unif();
  static double unif(double a, double b);
  static double norm();
  static double norm(double mean, double stddev);
  static double gamma(double k, double theta);

protected:
  static RNG *get();
  static RNG *_instance;
};

}

#endif // __SBB_RNG_HH__

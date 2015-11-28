#include "lib/randomvariable.hh"
#include "lib/operators.hh"
#include <iostream>


using namespace sbb;

int main(int argc, char *argv[]) {
  // stimulus at t=0:
  RandomVariable stimulus = GenericRandomVariable::delta(0);
  // add a 100 ms delay
  RandomVariable delay    = GenericRandomVariable::delta(100);
  // process stimulus
  RandomVariable process  = GenericRandomVariable::norm(250, 100);
  // some noise (+/- 40ms)
  RandomVariable noise    = GenericRandomVariable::norm(0, 20);

  // Assemble process
  RandomVariable out = (stimulus + delay + process + noise);

  Eigen::VectorXd pdf(1000);
  out.density().eval(0, 500, pdf);
  std::cout << pdf;

  return 0;
}

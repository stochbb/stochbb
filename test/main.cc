#include "lib/unittest.hh"
#include <iostream>

#include "convolutiontest.hh"
#include "compoundtest.hh"
#include "exactsamplertest.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;


int main(int argc, char *argv[]) {

  // Construct test-runner
  TestRunner runner(std::cout);

  // Add suites
  runner.addSuite(ConvolutionTest::suite());
  runner.addSuite(CompoundTest::suite());
  runner.addSuite(ExactSamplerTest::suite());

  // Exec tests:
  runner();

  return 0;
}

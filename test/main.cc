#include "lib/unittest.hh"
#include <iostream>

#include "convolutiontest.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;


int main(int argc, char *argv[]) {

  // Construct test-runner
  TestRunner runner(std::cout);

  // Add suites
  runner.addSuite(ConvolutionTest::suite());

  // Exec tests:
  runner();

  return 0;
}

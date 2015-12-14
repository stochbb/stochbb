#include "lib/unittest.hh"
#include <iostream>

using namespace stochbb::UnitTest;


int main(int argc, char *argv[]) {

  // Construct test-runner
  TestRunner runner(std::cout);

  // Add suites
  //runner.addSuite(...);

  // Exec tests:
  runner();

  return 0;
}

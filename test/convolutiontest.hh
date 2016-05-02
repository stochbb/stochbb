#ifndef CONVOLUTIONTEST_HH
#define CONVOLUTIONTEST_HH

#include "lib/unittest.hh"

namespace stochbb {

class ConvolutionTest : public UnitTest::TestCase
{
public:
  ConvolutionTest();

  void testAlignment();

public:
  static UnitTest::TestSuite *suite();

};

}

#endif // CONVOLUTIONTEST_HH

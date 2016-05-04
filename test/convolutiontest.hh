#ifndef CONVOLUTIONTEST_HH
#define CONVOLUTIONTEST_HH

#include "lib/unittest.hh"

namespace stochbb {

class ConvolutionTest : public UnitTest::TestCase
{
public:
  ConvolutionTest();

  void testAlignment();
  void testDeltaReduction();
  void testNormalReduction();
  void testGammaReduction();

public:
  static UnitTest::TestSuite *suite();

};

}

#endif // CONVOLUTIONTEST_HH

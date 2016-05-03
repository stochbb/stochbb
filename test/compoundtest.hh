#ifndef COMPOUNDTEST_HH
#define COMPOUNDTEST_HH

#include "lib/unittest.hh"

namespace stochbb {

class CompoundTest: public UnitTest::TestCase
{
public:
  CompoundTest();

  void testNormalCompound();
  void testGammaCompound();

public:
  static UnitTest::TestSuite *suite();
};

}

#endif // COMPOUNDTEST_HH

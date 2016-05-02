#ifndef COMPOUNDTEST_HH
#define COMPOUNDTEST_HH

#include "lib/unittest.hh"

namespace stochbb {

class CompoundTest: public UnitTest::TestCase
{
public:
  CompoundTest();

  void testNormalCompound();

public:
  static UnitTest::TestSuite *suite();
};

}

#endif // COMPOUNDTEST_HH

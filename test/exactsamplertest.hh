#ifndef EXACTSAMPLERTEST_HH
#define EXACTSAMPLERTEST_HH

#include "lib/unittest.hh"


namespace stochbb {

class ExactSamplerTest: public UnitTest::TestCase
{
public:
  ExactSamplerTest();

  void testChain();
  void testGammaCompound();

public:
  static UnitTest::TestSuite *suite();
};

}

#endif // EXACTSAMPLERTEST_HH

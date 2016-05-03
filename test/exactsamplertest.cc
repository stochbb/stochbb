#include "exactsamplertest.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;


ExactSamplerTest::ExactSamplerTest()
{
  // pass...
}

void
ExactSamplerTest::testGammaCompound() {

}


TestSuite *
ExactSamplerTest::suite() {
  TestSuite *suite = new TestSuite("ExactSampler");
  suite->addTest(new TestCaller<ExactSamplerTest>("gammaCompound", &ExactSamplerTest::testGammaCompound));
  return suite;
}





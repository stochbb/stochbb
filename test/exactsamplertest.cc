#include "exactsamplertest.hh"
#include "lib/exactsampler.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;


ExactSamplerTest::ExactSamplerTest()
{
  // pass...
}

void
ExactSamplerTest::testGammaCompound() {
  Var k = uniform(0,4);
  Var X = gamma(5*k+5, 10);

  size_t N = 100;
  Eigen::VectorXd dX(N);

  ExactSampler sampler(X);
  sampler.sample(dX);
}


TestSuite *
ExactSamplerTest::suite() {
  TestSuite *suite = new TestSuite("ExactSampler");
  suite->addTest(new TestCaller<ExactSamplerTest>("gammaCompound", &ExactSamplerTest::testGammaCompound));
  return suite;
}





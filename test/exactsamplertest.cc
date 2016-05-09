#include "exactsamplertest.hh"
#include "lib/exactsampler.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;


ExactSamplerTest::ExactSamplerTest()
{
  // pass...
}

void
ExactSamplerTest::testChain() {
  Var X = normal(0,1);
  Var Y = normal(0,1);
  Var Z = X + Y;
  size_t N = 100;
  Eigen::VectorXd dX(N);
  ExactSampler sampler(X);
  sampler.sample(dX);
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
  suite->addTest(new TestCaller<ExactSamplerTest>("normal chain sampler", &ExactSamplerTest::testChain));
  suite->addTest(new TestCaller<ExactSamplerTest>("gamma compound sampler", &ExactSamplerTest::testGammaCompound));
  return suite;
}





SET(STOCHBB_TEST_SOURCES main.cc convolutiontest.cc compoundtest.cc)

add_executable(stochbb_test ${STOCHBB_TEST_SOURCES})
target_link_libraries(stochbb_test ${LIBS} libsbb)

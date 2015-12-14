SET(STOCHBB_TEST_SOURCES main.cc )

add_executable(stochbb_test ${STOCHBB_TEST_SOURCES})
target_link_libraries(stochbb_test ${LIBS} libsbb)

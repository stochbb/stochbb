SET(stochbb_SOURCES main.cc)
SET(stochbb_HEADERS )

add_executable(stochbb ${stochbb_SOURCES})
target_link_libraries(stochbb ${LIBS} sbb)

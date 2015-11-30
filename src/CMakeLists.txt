SET(stochbb_SOURCES main.cc qcustomplot.cc plot.cc)
SET(stochbb_HEADERS plot.hh)
SET(stochbb_MOC_HEADERS qcustomplot.hh)

qt5_wrap_cpp(stochbb_MOC_SOURCES ${stochbb_MOC_HEADERS})
add_executable(stochbb ${stochbb_SOURCES} ${stochbb_MOC_SOURCES})
target_link_libraries(stochbb ${LIBS} sbb)

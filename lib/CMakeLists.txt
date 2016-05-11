SET(libsbb_SOURCES
    object.cc density.cc randomvariable.cc rng.cc chain.cc minmax.cc distribution.cc
    exception.cc math.cc api.cc option_parser.cc logger.cc operators.cc marginalsampler.cc
    exactsampler.cc mixture.cc compound.cc affinetrafo.cc conditional.cc unittest.cc cputime.cc
    reduction.cc)
SET(libsbb_HEADERS stochbb.hh
    object.hh density.hh randomvariable.hh rng.hh chain.hh minmax.hh distribution.hh
    exception.hh math.hh api.hh option_parser.hh logger.hh operators.hh marginalsampler.hh
    exactsampler.hh mixture.hh compound.hh affinetrafo.hh conditional.hh unittest.hh cputime.hh
    reduction.hh)

add_custom_target(sbb_headers SOURCES ${libsbb_HEADERS})

add_library(libsbb SHARED ${libsbb_SOURCES})
set_target_properties(libsbb PROPERTIES OUTPUT_NAME stochbb)
set_target_properties(libsbb PROPERTIES DEPENDS libsbb_HEADERS)
set_target_properties(libsbb PROPERTIES VERSION
 "${libsbb_VERSION_MAJOR}.${libsbb_VERSION_MINOR}.${libsbb_VERSION_PATCH}")
set_target_properties(libsbb PROPERTIES SOVERION ${libsbb_VERSION_MAJOR})
set_target_properties(libsbb PROPERTIES MACOSX_RPATH ${CMAKE_INSTALL_RPATH})
set_target_properties(libsbb PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_RPATH})
target_link_libraries(libsbb ${LIBS})

install(TARGETS libsbb DESTINATION "${CMAKE_INSTALL_LIBDIR}")
INSTALL_HEADERS_WITH_DIRECTORY("${libsbb_HEADERS}" "${CMAKE_INSTALL_INCLUDEDIR}/stochbb")

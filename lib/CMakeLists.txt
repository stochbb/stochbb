SET(libsbb_SOURCES
    object.cc density.cc randomvariable.cc rng.cc chain.cc minmax.cc simulation.cc xmlparser.cc
    exception.cc math.cc)
SET(libsbb_HEADERS mersennetwister.hh operators.hh
    object.hh density.hh randomvariable.hh rng.hh chain.hh minmax.hh simulation.hh xmlparser.hh
    exception.hh math.hh)

add_custom_target(sbb_headers SOURCES ${libsbb_HEADERS})

add_library(sbb SHARED ${libsbb_SOURCES})
target_link_libraries(sbb ${LIBS})

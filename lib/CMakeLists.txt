SET(libsbb_SOURCES
    object.cc density.cc randomvariable.cc rng.cc chain.cc minmax.cc simulation.cc xmlparser.cc
    exception.cc math.cc api.cc option_parser.cc logger.cc operators.cc marginalsampler.cc
    exactsampler.cc mixture.cc)
SET(libsbb_HEADERS mersennetwister.hh stochbb.hh
    object.hh density.hh randomvariable.hh rng.hh chain.hh minmax.hh simulation.hh xmlparser.hh
    exception.hh math.hh api.hh option_parser.hh logger.hh operators.hh marginalsampler.hh
    exactsampler.hh mixture.hh)

add_custom_target(sbb_headers SOURCES ${libsbb_HEADERS})

add_library(sbb SHARED ${libsbb_SOURCES})
target_link_libraries(sbb ${LIBS})

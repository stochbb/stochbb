SET(example_FILES simple_chain.xml gamma.xml minmax.xml exp.xml exgauss.xml chisq.xml mixture.xml
    compound.xml example.xml)
add_custom_target(examples SOURCES ${example_FILES})

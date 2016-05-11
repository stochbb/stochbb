SET(example_FILES condchain.py condchain.R chain.py chain.R pod.py pod.R
    range_estimate.py range_estimate.R directConvolution.py interval.py)
add_custom_target(examples SOURCES ${example_FILES})

SET(R_module_SOURCES src/module.cc)
SET(R_module_SCRIPTS R/stochbb.R)
SET(R_module_DOC_FILES man/stochbb-package.Rd man/delta.Rd man/uniform.Rd man/normal.Rd
    man/gamma.Rd man/invgamma.Rd man/weibull.Rd man/chain.Rd man/affine.Rd man/minimum.Rd
    man/maximum.Rd man/independent.Rd man/mixture.Rd man/conditional.Rd man/condchain.Rd
    man/atomic.Rd man/derived.Rd man/studt.Rd
    man/Container-class.Rd man/Var-class.Rd man/Density-class.Rd man/ExactSampler-class.Rd)
SET(R_module_FILES ${R_module_SOURCES} ${R_module_SCRIPTS} ${R_module_DOC_FILES})
add_custom_target(R_module SOURCES ${R_module_FILES})

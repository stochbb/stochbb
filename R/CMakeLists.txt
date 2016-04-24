SET(R_module_SOURCES src/module.cc)
SET(R_module_SCRIPTS R/stochbb.R)
SET(R_module_DOC_FILES man/stochbb-package.Rd man/Var.Rd man/uniform.Rd man/normal.Rd man/gamma.Rd
    man/invgamma.Rd man/weibull.Rd)
SET(R_module_FILES ${R_module_SOURCES} ${R_module_SCRIPTS} ${R_module_DOC_FILES})
add_custom_target(R_module SOURCES ${R_module_FILES})

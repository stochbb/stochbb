# - Find Sundials NVector
#
#  NVECTOR_INCLUDES    - where to find nvector/nvector_serial.h
#  NVECTOR_LIBRARIES   - List of libraries when using NVector.
#  NVECTOR_FOUND       - True if Sundials NVector was found.

if (NVECTOR_INCLUDES)
  # Already in cache, be silent
  set (NVECTOR_FIND_QUIETLY TRUE)
endif (NVECTOR_INCLUDES)

find_path(NVECTOR_INCLUDE_DIRS nvector/nvector_serial.h)
find_library(NVECTOR_LIBRARIES NAMES sundials_nvecserial)

# handle the QUIETLY and REQUIRED arguments and set NVECTOR_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NVECTOR DEFAULT_MSG NVECTOR_LIBRARIES NVECTOR_INCLUDE_DIRS)

mark_as_advanced (NVECTOR_LIBRARIES NVECTOR_INCLUDE_DIRS)

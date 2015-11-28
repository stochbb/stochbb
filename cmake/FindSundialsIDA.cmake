# - Find Sundials IDA solver
#
#  IDA_INCLUDES    - where to find ida/ida.h
#  IDA_LIBRARIES   - List of libraries when using IDA.
#  IDA_FOUND       - True if Sundials IDA was found.

if (IDA_INCLUDES)
  # Already in cache, be silent
  set (IDA_FIND_QUIETLY TRUE)
endif (IDA_INCLUDES)

find_path(IDA_INCLUDE_DIRS ida/ida.h)
find_library(IDA_LIBRARIES NAMES sundials_ida)

# handle the QUIETLY and REQUIRED arguments and set IDA_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (IDA DEFAULT_MSG IDA_LIBRARIES IDA_INCLUDE_DIRS)

mark_as_advanced (IDA_LIBRARIES IDA_INCLUDE_DIRS)

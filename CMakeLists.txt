# CMake script for step-85 project
# Based on bare-dealii-app template

# Set the name of the project and target
SET(_library_target step-85-lib)

# Declare all source files the target consists of:
file(GLOB _files source/*.cc include/*.h)
file(GLOB _execs execs/*.cc)

INCLUDE_DIRECTORIES(include)

CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

FIND_PACKAGE(deal.II 9.5 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(${_library_target})

# We add one library and one target for each type of deal.II library
# we found. If you compiled deal.II with both Release and Debug
# mode, this will generate both Release and Debug programs for you.
SET(_d2_build_types "Release;Debug")
SET(Release_postfix "")
SET(Debug_postfix ".g")

FOREACH(_build_type ${_d2_build_types})
  # Postfix to use everywhere
  SET(_p "${${_build_type}_postfix}")

  # Only build this type, if deal.II was compiled with it.
  IF(CMAKE_BUILD_TYPE MATCHES "${_build_type}" AND
    DEAL_II_BUILD_TYPE MATCHES "${_build_type}")
    MESSAGE("-- Found ${_build_type} version of deal.II.")

    STRING(TOUPPER "${_build_type}" _BUILD_TYPE)

    # Setup library
    SET(_lib "${_library_target}${_p}")
    MESSAGE("-- Configuring library ${_lib}")

    ADD_LIBRARY(${_lib} SHARED ${_files})
    DEAL_II_SETUP_TARGET(${_lib} ${_BUILD_TYPE})

    # Setup executables
    FOREACH(__exec ${_execs})
      get_filename_component(exec_name ${__exec} NAME_WE)
      SET(_exe "${exec_name}${_p}")
      MESSAGE("-- Configuring executable ${_exe}")
      ADD_EXECUTABLE(${_exe} ${__exec})
      TARGET_LINK_LIBRARIES(${_exe} ${_lib})
      DEAL_II_SETUP_TARGET(${_exe} ${_BUILD_TYPE})
    ENDFOREACH()
  ENDIF()
ENDFOREACH()

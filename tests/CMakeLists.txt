INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

# For 2D tests, use 2D library
SET(TEST_LIBRARIES_DEBUG step-85_2d.g)
SET(TEST_LIBRARIES_RELEASE step-85_2d)

# Note: DEAL_II_PICKUP_TESTS() picks up all *.cc files in the directory
# For now, we manually handle 3D tests below
# List of 2D tests (those that should link against 2D library)
SET(_2d_tests
  boundary_passes
  full_residual_patches
  full_residual_patches_passes
  shy_patches
)

# List of 3D tests (those that should link against 3D library)
SET(_3d_tests
  full_residual_patches_3D
  shy_patches_3D
)

# Build 2D tests
FOREACH(_test ${_2d_tests})
  SET(_test_source "${CMAKE_CURRENT_SOURCE_DIR}/${_test}.cc")
  IF(EXISTS ${_test_source})
    STRING(TOUPPER "${CMAKE_BUILD_TYPE}" _build)
    IF(_build MATCHES "DEBUG")
      ADD_EXECUTABLE(${_test}.g ${_test_source})
      DEAL_II_SETUP_TARGET(${_test}.g DEBUG)
      TARGET_LINK_LIBRARIES(${_test}.g ${TEST_LIBRARIES_DEBUG})
    ENDIF()
    IF(_build MATCHES "RELEASE")
      ADD_EXECUTABLE(${_test} ${_test_source})
      DEAL_II_SETUP_TARGET(${_test} RELEASE)
      TARGET_LINK_LIBRARIES(${_test} ${TEST_LIBRARIES_RELEASE})
    ENDIF()
  ENDIF()
ENDFOREACH()

# Build 3D tests with 3D library
SET(TEST_LIBRARIES_DEBUG step-85_3d.g)
SET(TEST_LIBRARIES_RELEASE step-85_3d)

FOREACH(_test ${_3d_tests})
  SET(_test_source "${CMAKE_CURRENT_SOURCE_DIR}/${_test}.cc")
  IF(EXISTS ${_test_source})
    STRING(TOUPPER "${CMAKE_BUILD_TYPE}" _build)
    IF(_build MATCHES "DEBUG")
      ADD_EXECUTABLE(${_test}.g ${_test_source})
      DEAL_II_SETUP_TARGET(${_test}.g DEBUG)
      TARGET_LINK_LIBRARIES(${_test}.g ${TEST_LIBRARIES_DEBUG})
    ENDIF()
    IF(_build MATCHES "RELEASE")
      ADD_EXECUTABLE(${_test} ${_test_source})
      DEAL_II_SETUP_TARGET(${_test} RELEASE)
      TARGET_LINK_LIBRARIES(${_test} ${TEST_LIBRARIES_RELEASE})
    ENDIF()
  ENDIF()
ENDFOREACH()

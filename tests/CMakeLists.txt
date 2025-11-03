INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

# Link tests against the 2D library (tests are currently 2D-only)
SET(TEST_LIBRARIES_DEBUG step-85_2d.g)
SET(TEST_LIBRARIES_RELEASE step-85_2d)

DEAL_II_PICKUP_TESTS()

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

# Link tests against the library
SET(TEST_LIBRARIES_DEBUG step-85-lib.g)

DEAL_II_PICKUP_TESTS()

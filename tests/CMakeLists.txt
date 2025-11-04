INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

# Setup for 2D tests
SET(TEST_LIBRARIES_DEBUG step-85_2d.g)
SET(TEST_LIBRARIES_RELEASE step-85_2d)

# Only pick up 2D tests (explicitly list them) - need to include "tests/" prefix
SET(TEST_PICKUP_REGEX "^tests/(boundary_passes|full_residual_patches|full_residual_patches_passes|shy_patches)$")

# Pickup all 2D tests
DEAL_II_PICKUP_TESTS()

# Now handle 3D tests separately
# Only pick up 3D tests (full_residual_patches_3D only - shy_patches don't work in 3D)
SET(TEST_PICKUP_REGEX "^tests/(full_residual_patches_3D)$")
SET(TEST_LIBRARIES_DEBUG step-85_3d.g)
SET(TEST_LIBRARIES_RELEASE step-85_3d)

# Pickup 3D tests
DEAL_II_PICKUP_TESTS()

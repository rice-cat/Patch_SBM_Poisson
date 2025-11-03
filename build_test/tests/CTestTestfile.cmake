# CMake generated Testfile for 
# Source directory: /workspace/tests
# Build directory: /workspace/build_test/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(tests/full_residual_patches.debug "/usr/bin/cmake" "-DTRGT=tests.full_residual_patches.debug.test" "-DTEST=tests/full_residual_patches.debug" "-DEXPECT=PASSED" "-DBINARY_DIR=/workspace/build_test" "-P" "/usr/local/share/deal.II/scripts/run_test.cmake")
set_tests_properties(tests/full_residual_patches.debug PROPERTIES  LABEL "tests" PROCESSORS "1" TIMEOUT "600" WORKING_DIRECTORY "/workspace/build_test/tests/full_residual_patches.debug/serial" _BACKTRACE_TRIPLES "/usr/local/share/deal.II/macros/macro_deal_ii_add_test.cmake;621;add_test;/usr/local/share/deal.II/macros/macro_deal_ii_pickup_tests.cmake;383;deal_ii_add_test;/workspace/tests/CMakeLists.txt;7;DEAL_II_PICKUP_TESTS;/workspace/tests/CMakeLists.txt;0;")
add_test(tests/shy_patches.debug "/usr/bin/cmake" "-DTRGT=tests.shy_patches.debug.test" "-DTEST=tests/shy_patches.debug" "-DEXPECT=PASSED" "-DBINARY_DIR=/workspace/build_test" "-P" "/usr/local/share/deal.II/scripts/run_test.cmake")
set_tests_properties(tests/shy_patches.debug PROPERTIES  LABEL "tests" PROCESSORS "1" TIMEOUT "600" WORKING_DIRECTORY "/workspace/build_test/tests/shy_patches.debug/serial" _BACKTRACE_TRIPLES "/usr/local/share/deal.II/macros/macro_deal_ii_add_test.cmake;621;add_test;/usr/local/share/deal.II/macros/macro_deal_ii_pickup_tests.cmake;383;deal_ii_add_test;/workspace/tests/CMakeLists.txt;7;DEAL_II_PICKUP_TESTS;/workspace/tests/CMakeLists.txt;0;")

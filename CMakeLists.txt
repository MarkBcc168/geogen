cmake_minimum_required(VERSION 3.16)
project(geogen VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(geogen src/geogen.cpp)
target_compile_options(geogen PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic -Wconversion>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>)

include(CTest)
if(BUILD_TESTING)
  add_test(NAME sample_generate
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/midpoint_perpendicular.geogen)
  set_tests_properties(sample_generate PROPERTIES
    PASS_REGULAR_EXPRESSION "GEOGEN REPORT")
  add_test(NAME sample_prove
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/orthocenter.geogen)
  set_tests_properties(sample_prove PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED perpendicular")
  add_test(NAME sample_known_intersection
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/cyclic_second_intersection.geogen)
  set_tests_properties(sample_known_intersection PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic")
  add_test(NAME sample_unproved_conjecture
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/unproved_concyclic.geogen)
  set_tests_properties(sample_unproved_conjecture PROPERTIES
    PASS_REGULAR_EXPRESSION "NONTRIVIAL concyclic")
  add_test(NAME no_illegal_angle_halving
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/no_angle_halving.geogen)
  set_tests_properties(no_illegal_angle_halving PROPERTIES
    PASS_REGULAR_EXPRESSION "UNPROVED parallel")
  add_test(NAME adjustable_automatic_points
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/automatic_points.geogen)
  set_tests_properties(adjustable_automatic_points PROPERTIES
    PASS_REGULAR_EXPRESSION "points=10")
  add_test(NAME discard_coincident_points
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/coincident_points.geogen)
  set_tests_properties(discard_coincident_points PROPERTIES
    PASS_REGULAR_EXPRESSION "points=5"
    FAIL_REGULAR_EXPRESSION "POINT X ")
endif()

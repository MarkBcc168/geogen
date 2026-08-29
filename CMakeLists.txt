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
    PASS_REGULAR_EXPRESSION "points=50")
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
    PASS_REGULAR_EXPRESSION "stable_coincidences=0")
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
  add_test(NAME right_triangle_midpoint_rule
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/right_triangle_midpoint.geogen)
  set_tests_properties(right_triangle_midpoint_rule PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED equal_distance.M,A,M,D")
  add_test(NAME prompt_concyclic_filter
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/prompt_concyclic.geogen)
  set_tests_properties(prompt_concyclic_filter PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.A,E,F,H")
  add_test(NAME definition_collinearity_filter
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/definition_collinearity.geogen)
  set_tests_properties(definition_collinearity_filter PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED collinear")
  add_test(NAME midpoint_midline_rule
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/midline.geogen)
  set_tests_properties(midpoint_midline_rule PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED parallel")
  add_test(NAME automatic_known_root_line_circle
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/automatic_line_circle.geogen)
  set_tests_properties(automatic_known_root_line_circle PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT X.secant,omega,A")
  add_test(NAME circumcenter_isosceles_filter
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/circumcenter_isosceles.geogen)
  set_tests_properties(circumcenter_isosceles_filter PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.A,O1,O2,O3")
  add_test(NAME perpendicular_bisector_circumcenter_filter
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/circumcenter_isosceles.geogen)
  set_tests_properties(perpendicular_bisector_circumcenter_filter PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.B,F,M,O1")
  add_test(NAME center_locus_collinearity_filter
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/center_loci.geogen)
  set_tests_properties(center_locus_collinearity_filter PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED collinear"
    FAIL_REGULAR_EXPRESSION "UNPROVED collinear")
  add_test(NAME ancestry_scoped_filter
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/ancestry_scope.geogen)
  set_tests_properties(ancestry_scoped_filter PROPERTIES
    PASS_REGULAR_EXPRESSION "summary stable_coincidences=")
  add_test(NAME random_generation_reaches_cap
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/random_generation.geogen)
  set_tests_properties(random_generation_reaches_cap PROPERTIES
    PASS_REGULAR_EXPRESSION "points=100")
  add_test(NAME random_generation_line_line
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/random_generation.geogen)
  set_tests_properties(random_generation_line_line PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT Xll.* = intersect")
  add_test(NAME random_generation_known_root_circles
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/random_generation.geogen)
  set_tests_properties(random_generation_known_root_circles PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT X(lc|cc).* = other_intersection")
  add_test(NAME point_output_definitions
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/output_definitions.geogen)
  set_tests_properties(point_output_definitions PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT X = reflect.A,line.B,C")
  add_test(NAME input_points_stay_atomic
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/output_definitions.geogen)
  set_tests_properties(input_points_stay_atomic PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT Y = reflect.C,E")
  add_test(NAME orthocenter_reflection_cyclic
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/orthocenter_reflection.geogen)
  set_tests_properties(orthocenter_reflection_cyclic PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.A,B,C,D")
  add_test(NAME orthocenter_reflection_collinear
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/orthocenter_reflection.geogen)
  set_tests_properties(orthocenter_reflection_collinear PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED collinear.A,H,D")
  add_test(NAME shared_midpoint_parallelogram
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/parallelogram_reflection.geogen)
  set_tests_properties(shared_midpoint_parallelogram PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED parallel.A,B,C,D")
  add_test(NAME point_reflection_midpoint_cyclic
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/parallelogram_reflection.geogen)
  set_tests_properties(point_reflection_midpoint_cyclic PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.B,H,C,D")
  add_test(NAME definition_line_incidences
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/definition_incidences.geogen)
  set_tests_properties(definition_line_incidences PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED collinear.A,M,D")
  add_test(NAME definition_circle_incidences
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/definition_incidences.geogen)
  set_tests_properties(definition_circle_incidences PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.A,B,C,D")
  add_test(NAME incenter_foot_equal_tangents
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/incenter_kite.geogen)
  set_tests_properties(incenter_foot_equal_tangents PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED equal_distance.A,E,A,F")
  add_test(NAME incenter_foot_kite
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/incenter_kite.geogen)
  set_tests_properties(incenter_foot_kite PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED perpendicular.A,I,E,F")
  add_test(NAME cyclic_radius_completion
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/cyclic_radius_completion.geogen)
  set_tests_properties(cyclic_radius_completion PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED equal_distance.O,A,O,D")
  add_test(NAME maximal_definition_collinearity
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/maximal_definition_collinearity.geogen)
  set_tests_properties(maximal_definition_collinearity PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED collinear.A,B,F,P,Q,R,S,W")
  add_test(NAME circumcenter_chord_midpoint
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/circumcenter_chord_midpoint.geogen)
  set_tests_properties(circumcenter_chord_midpoint PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED equal_distance.A,M,M,B")
endif()

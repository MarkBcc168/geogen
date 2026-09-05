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
    PASS_REGULAR_EXPRESSION "generated_points=314")
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
    PASS_REGULAR_EXPRESSION "combined_coincidences=0")
  add_test(NAME no_illegal_angle_halving
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/no_angle_halving.geogen)
  set_tests_properties(no_illegal_angle_halving PROPERTIES
    PASS_REGULAR_EXPRESSION "UNPROVED parallel")
  add_test(NAME adjustable_automatic_points
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/automatic_points.geogen)
  set_tests_properties(adjustable_automatic_points PROPERTIES
    PASS_REGULAR_EXPRESSION "generated_points=24")
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
    PASS_REGULAR_EXPRESSION "generated_points=5")
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
    PASS_REGULAR_EXPRESSION "summary combined_coincidences=")
  add_test(NAME random_generation_reaches_cap
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/random_generation.geogen)
  set_tests_properties(random_generation_reaches_cap PROPERTIES
    PASS_REGULAR_EXPRESSION "generated_points=100")
  add_test(NAME random_generation_line_line
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/random_generation.geogen)
  set_tests_properties(random_generation_line_line PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT Xll.* = intersect")
  add_test(NAME random_generation_known_root_circles
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/random_generation.geogen)
  set_tests_properties(random_generation_known_root_circles PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT X(lc|cc).* = other_intersection")
  add_test(NAME random_generation_uses_circumcircles
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/random_generation.geogen)
  set_tests_properties(random_generation_uses_circumcircles PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT X(lc|cc).*circumcircle")
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
  add_test(NAME iran_lemma_diameter_circle
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/iran_lemma.geogen)
  set_tests_properties(iran_lemma_diameter_circle PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.B,I,X,D,F")
  add_test(NAME iran_lemma_right_angle
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/iran_lemma.geogen)
  set_tests_properties(iran_lemma_right_angle PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED perpendicular.A,X,B,X")
  add_test(NAME isosceles_incenter_symmetry
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/stress_isosceles_incenter.geogen)
  set_tests_properties(isosceles_incenter_symmetry PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED collinear.A,B,H")
  add_test(NAME angle_bisector_reflection_symmetry
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/stress_bisector_reflection.geogen)
  set_tests_properties(angle_bisector_reflection_symmetry PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED collinear.A,C,Q,M")
  add_test(NAME circle_circle_center_axis
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/stress_circle_circle_axis.geogen)
  set_tests_properties(circle_circle_center_axis PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED perpendicular.A,X,B,C")
  add_test(NAME delayed_reflection_axis_incidence
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/stress_reflected_circumcircle.geogen)
  set_tests_properties(delayed_reflection_axis_incidence PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.A,B,C,D")
  add_test(NAME affine_nested_reflections
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/affine_nested_reflections.geogen)
  set_tests_properties(affine_nested_reflections PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED collinear.C,R,M")
  add_test(NAME affine_implicit_midpoint
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/affine_nested_reflections.geogen)
  set_tests_properties(affine_implicit_midpoint PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED equal_distance.C,T,T,R")
  add_test(NAME affine_facts_feed_angle_chase
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/affine_nested_reflections.geogen)
  set_tests_properties(affine_facts_feed_angle_chase PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED parallel.T,N,R,A")
  add_test(NAME affine_line_intersection
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/affine_median_intersection.geogen)
  set_tests_properties(affine_line_intersection PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED collinear.A,K,G")
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
  add_test(NAME coincident_lines_are_canonical
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/coincident_lines.geogen)
  set_tests_properties(coincident_lines_are_canonical PROPERTIES
    PASS_REGULAR_EXPRESSION "lines=3")
  add_test(NAME coincident_line_output_uses_carrier
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/coincident_lines.geogen)
  set_tests_properties(coincident_line_output_uses_carrier PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT X = foot.A,line.B,C")
  add_test(NAME reflected_homothety_circumcenter_circle
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/circumcenter_reflection_homothety.geogen)
  set_tests_properties(reflected_homothety_circumcenter_circle PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.A,B,C,O")
  add_test(NAME incenter_side_circumcenter_circle
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/incenter_side_circumcenter.geogen)
  set_tests_properties(incenter_side_circumcenter_circle PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.A,B,C,O")
  add_test(NAME perpendicular_bisector_reflection_circle
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/perpendicular_bisector_reflection.geogen)
  set_tests_properties(perpendicular_bisector_reflection_circle PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.A,B,C,O1,R1,R2,O2")
  add_test(NAME orthogonal_isosceles_trapezoids
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/orthogonal_isosceles_trapezoids.geogen)
  set_tests_properties(orthogonal_isosceles_trapezoids PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.A,F,G,R2")
  add_test(NAME fixed_point_cross_prover_propagation
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/orthogonal_isosceles_trapezoids.geogen)
  set_tests_properties(fixed_point_cross_prover_propagation PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED equal_distance.O,A,O,R2")
  add_test(NAME circle_chord_reflection_rule
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/chord_reflection.geogen)
  set_tests_properties(circle_chord_reflection_rule PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.E,P,H,H1")
  add_test(NAME reflection_halfturn_composition_rule
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/reflection_halfturn.geogen)
  set_tests_properties(reflection_halfturn_composition_rule PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED concyclic.X,Y,C,D")
  add_test(NAME derived_orthocenter_third_altitude
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/derived_orthocenter.geogen)
  set_tests_properties(derived_orthocenter_third_altitude PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED perpendicular.H,X,A,Y")
  add_test(NAME formal_affine_unknown_point_parallelism
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/derived_orthocenter.geogen)
  set_tests_properties(formal_affine_unknown_point_parallelism PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED parallel.A,Y,X,P")
  add_test(NAME derived_orthocenter_propagates_to_parallel_side
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/derived_orthocenter.geogen)
  set_tests_properties(derived_orthocenter_propagates_to_parallel_side PROPERTIES
    PASS_REGULAR_EXPRESSION "PROVED perpendicular.H,X,X,P")
  add_test(NAME converse_cyclic_names_quadrilateral
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/circumcenter_isosceles.geogen)
  set_tests_properties(converse_cyclic_names_quadrilateral PROPERTIES
    PASS_REGULAR_EXPRESSION "converse cyclic angle theorem B,M,F,O1")
  add_test(NAME report_marks_unused_initial_points
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/unused_initial.geogen)
  set_tests_properties(report_marks_unused_initial_points PROPERTIES
    PASS_REGULAR_EXPRESSION "unused_initial=C")
  add_test(NAME report_marks_unused_explicit_input_points
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/unused_input_point.geogen)
  set_tests_properties(report_marks_unused_explicit_input_points PROPERTIES
    PASS_REGULAR_EXPRESSION "collinear.A,B,M.*unused_initial=C,E")
  add_test(NAME report_marks_isosceles_trapezoid
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/report_annotations.geogen)
  set_tests_properties(report_marks_isosceles_trapezoid PROPERTIES
    PASS_REGULAR_EXPRESSION "isosceles_trapezoid=")
  add_test(NAME report_marks_two_right_angles
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/report_annotations.geogen)
  set_tests_properties(report_marks_two_right_angles PROPERTIES
    PASS_REGULAR_EXPRESSION "two_right_angles=")
  add_test(NAME report_omits_auxiliary_point_definitions
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/report_annotations.geogen)
  set_tests_properties(report_omits_auxiliary_point_definitions PROPERTIES
    FAIL_REGULAR_EXPRESSION "POINT I =")
  add_test(NAME symmetric_generation_pairs_points
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/symmetric_generation.geogen)
  set_tests_properties(symmetric_generation_pairs_points PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT S\\$Foot\\$4 = foot.C,perpendicular_bisector.C,A")
  add_test(NAME symmetric_generation_marks_pairs
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/symmetric_generation.geogen)
  set_tests_properties(symmetric_generation_marks_pairs PROPERTIES
    PASS_REGULAR_EXPRESSION "symmetry=asymmetric.*symmetric_partner=collinear.A,C,S\\$Foot\\$4")
  add_test(NAME symmetric_foot_construction
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/symmetric_foot.geogen)
  set_tests_properties(symmetric_foot_construction PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT S\\$F = foot.C,line.A,B")
  add_test(NAME symmetric_foot_coincidence_pair
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/symmetric_foot.geogen)
  set_tests_properties(symmetric_foot_coincidence_pair PROPERTIES
    PASS_REGULAR_EXPRESSION "symmetry=asymmetric.*symmetric_partner=collinear.A,B,S\\$F")
  add_test(NAME show_all_generated_points
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/show_all_points.geogen)
  set_tests_properties(show_all_generated_points PROPERTIES
    PASS_REGULAR_EXPRESSION "generated_points=10 reported_points=10 hidden_points=0")
  add_test(NAME disable_point_reflection_generation
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/construction_controls.geogen)
  set_tests_properties(disable_point_reflection_generation PROPERTIES
    FAIL_REGULAR_EXPRESSION "POINT Rp\\$")
  add_test(NAME disable_incenter_generation
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/construction_controls.geogen)
  set_tests_properties(disable_incenter_generation PROPERTIES
    FAIL_REGULAR_EXPRESSION "POINT I\\$")
  add_test(NAME symmetric_only_filters_asymmetric_pairs
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/symmetric_only.geogen)
  set_tests_properties(symmetric_only_filters_asymmetric_pairs PROPERTIES
    PASS_REGULAR_EXPRESSION "EASY concyclic.B,C,F,S\\$F.*symmetry=self"
    FAIL_REGULAR_EXPRESSION "symmetry=asymmetric")
  add_test(NAME symmetric_only_filters_point_listing
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/symmetric_only.geogen)
  set_tests_properties(symmetric_only_filters_point_listing PROPERTIES
    PASS_REGULAR_EXPRESSION "generated_points=5 reported_points=4 hidden_points=1"
    FAIL_REGULAR_EXPRESSION "POINT A =")
  add_test(NAME independent_trials_combine_generated_graphs
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/combined_trials.geogen)
  set_tests_properties(independent_trials_combine_generated_graphs PROPERTIES
    PASS_REGULAR_EXPRESSION "generated_points=37 reported_points=37 hidden_points=0")
  add_test(NAME independent_trials_namespace_points
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/combined_trials.geogen)
  set_tests_properties(independent_trials_namespace_points PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT T1\\$")
  add_test(NAME independent_trials_second_namespace
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/combined_trials.geogen)
  set_tests_properties(independent_trials_second_namespace PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT T2\\$")
  add_test(NAME two_point_lines_generate_intersections
    COMMAND geogen ${CMAKE_CURRENT_SOURCE_DIR}/examples/two_point_line_generation.geogen)
  set_tests_properties(two_point_lines_generate_intersections PROPERTIES
    PASS_REGULAR_EXPRESSION "POINT Xll\\$.* = intersect.line.C,B.,line.D,A")
endif()

cmake_minimum_required(VERSION 3.20)
project(femto_same_mixed_comparison LANGUAGES CXX)

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
add_executable(compare_same_mixed main.cpp Comparison.cpp Plotting.cpp)
target_compile_features(compare_same_mixed PRIVATE cxx_std_17)
target_compile_options(compare_same_mixed PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(compare_same_mixed PRIVATE ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad)
set_target_properties(compare_same_mixed PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

enable_testing()
add_test(NAME comparison_math_and_flow_contract COMMAND compare_same_mixed --self-test)

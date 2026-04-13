cmake_minimum_required(VERSION 3.16)

project(EventgenFemto3D LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

find_package(ROOT REQUIRED COMPONENTS Graf Hist MathCore Tree TreePlayer RIO Minuit2)
find_package(tomlplusplus CONFIG REQUIRED
  HINTS
    /opt/homebrew/lib/cmake/tomlplusplus
    /opt/homebrew/opt/tomlplusplus
)

add_library(eventgen_femto_3d_core
  src/Config.cpp
  src/EventPlane.cpp
  src/Histogramming.cpp
  src/InputReader.cpp
  src/ProjectionFit.cpp
  src/Source_extraction.cpp
  src/Workflow.cpp
)

target_include_directories(eventgen_femto_3d_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(eventgen_femto_3d_core
  PUBLIC
    ROOT::Graf
    ROOT::Hist
    ROOT::MathCore
    ROOT::Minuit2
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::RIO
    tomlplusplus::tomlplusplus
)

add_executable(eventgen_femto_3d
  src/main.cpp
)

target_link_libraries(eventgen_femto_3d
  PRIVATE
    eventgen_femto_3d_core
)

add_executable(directional_covariance_failure_semantics_test
  tests/DirectionalCovarianceFailureSemanticsTest.cpp
)

target_include_directories(directional_covariance_failure_semantics_test
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(directional_covariance_failure_semantics_test
  PRIVATE
    ROOT::Graf
    ROOT::Hist
    ROOT::MathCore
    ROOT::Minuit2
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::RIO
)

add_test(
  NAME directional_covariance_failure_semantics_test
  COMMAND directional_covariance_failure_semantics_test
)

add_executable(alpha_hbt_error_semantics_test
  tests/AlphaHbtErrorSemanticsTest.cpp
)

target_include_directories(alpha_hbt_error_semantics_test
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(alpha_hbt_error_semantics_test
  PRIVATE
    ROOT::Graf
    ROOT::Hist
    ROOT::MathCore
    ROOT::Minuit2
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::RIO
)

add_test(
  NAME alpha_hbt_error_semantics_test
  COMMAND alpha_hbt_error_semantics_test
)

add_executable(fit_summary_directional_error_mask_test
  tests/FitSummaryDirectionalErrorMaskTest.cpp
)

target_include_directories(fit_summary_directional_error_mask_test
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(fit_summary_directional_error_mask_test
  PRIVATE
    ROOT::Graf
    ROOT::Hist
    ROOT::MathCore
    ROOT::Minuit2
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::RIO
)

add_test(
  NAME fit_summary_directional_error_mask_test
  COMMAND fit_summary_directional_error_mask_test
)

add_executable(r2_summary_policy_test
  tests/R2SummaryPolicyTest.cpp
)

target_link_libraries(r2_summary_policy_test
  PRIVATE
    eventgen_femto_3d_core
)

add_test(
  NAME r2_summary_policy_test
  COMMAND r2_summary_policy_test
)

add_executable(config_parse_validation_test
  tests/ConfigParseValidationTest.cpp
)

target_link_libraries(config_parse_validation_test
  PRIVATE
    eventgen_femto_3d_core
)

add_test(
  NAME config_parse_validation_test
  COMMAND config_parse_validation_test
)

add_executable(root_runtime_probe
  tests/RootRuntimeProbe.cpp
)

target_link_libraries(root_runtime_probe
  PRIVATE
    eventgen_femto_3d_core
)

add_executable(input_adapter_test
  tests/InputAdapterTest.cpp
)

target_link_libraries(input_adapter_test
  PRIVATE
    eventgen_femto_3d_core
)

add_test(
  NAME input_adapter_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe>
          $<TARGET_FILE:input_adapter_test>
)

add_executable(legacy_workflow_smoke_test
  tests/LegacyWorkflowSmokeTest.cpp
)

target_link_libraries(legacy_workflow_smoke_test
  PRIVATE
    eventgen_femto_3d_core
)

add_test(
  NAME legacy_workflow_smoke_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe>
          $<TARGET_FILE:legacy_workflow_smoke_test>
)

add_executable(r2_summary_visibility_test
  tests/R2SummaryVisibilityTest.cpp
)

target_link_libraries(r2_summary_visibility_test
  PRIVATE
    eventgen_femto_3d_core
)

add_test(
  NAME r2_summary_visibility_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe>
          $<TARGET_FILE:r2_summary_visibility_test>
)

set_tests_properties(
  input_adapter_test
  legacy_workflow_smoke_test
  r2_summary_visibility_test
  PROPERTIES
    SKIP_RETURN_CODE 77
)

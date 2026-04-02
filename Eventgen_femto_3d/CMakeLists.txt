cmake_minimum_required(VERSION 3.16)

project(EventgenFemto3D LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

find_package(ROOT REQUIRED COMPONENTS Graf Hist MathCore Tree TreePlayer RIO Minuit2)

add_executable(eventgen_femto_3d
  src/EventPlane.cpp
  src/Histogramming.cpp
  src/ProjectionFit.cpp
  src/Source_extraction.cpp
  src/Workflow.cpp
  src/main.cpp
)

target_include_directories(eventgen_femto_3d
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(eventgen_femto_3d
  PRIVATE
    ROOT::Graf
    ROOT::Hist
    ROOT::MathCore
    ROOT::Minuit2
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::RIO
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

cmake_minimum_required(VERSION 3.16)

project(ExpFemto3D LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json" FORCE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE PATH "Executable output directory")

if(CMAKE_CONFIGURATION_TYPES)
  foreach(config ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER "${config}" config_upper)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${config_upper}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endforeach()
endif()

find_package(Threads REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad Tree TreePlayer MathCore Minuit)
find_package(tomlplusplus CONFIG REQUIRED
  HINTS
    /opt/homebrew/lib/cmake/tomlplusplus
    /opt/homebrew/opt/tomlplusplus
)

option(EXP_FEMTO_3D_ENABLE_CATS "Enable optional CATS finite-source Coulomb support" ON)
set(EXP_FEMTO_3D_CATS_ROOT "/Users/allenzhou/Research_software/CATS" CACHE PATH "Local CATS installation root")
set(EXP_FEMTO_3D_GSL_ROOT "/opt/homebrew/opt/gsl" CACHE PATH "Local GSL installation root")
set(EXP_FEMTO_3D_HAS_CATS OFF)
if(EXP_FEMTO_3D_ENABLE_CATS)
  set(EXP_FEMTO_3D_CATS_INCLUDE_DIR "${EXP_FEMTO_3D_CATS_ROOT}/install/include")
  set(EXP_FEMTO_3D_CATS_LIBRARY "${EXP_FEMTO_3D_CATS_ROOT}/install/lib/libCATS.dylib")
  find_package(GSL QUIET HINTS "${EXP_FEMTO_3D_GSL_ROOT}")
  set(EXP_FEMTO_3D_GSL_INCLUDE_DIR "${EXP_FEMTO_3D_GSL_ROOT}/include")
  set(EXP_FEMTO_3D_GSL_LIBRARY "${EXP_FEMTO_3D_GSL_ROOT}/lib/libgsl.dylib")
  set(EXP_FEMTO_3D_GSLCBLAS_LIBRARY "${EXP_FEMTO_3D_GSL_ROOT}/lib/libgslcblas.dylib")
  set(EXP_FEMTO_3D_GSL_AVAILABLE ${GSL_FOUND})
  if(NOT EXP_FEMTO_3D_GSL_AVAILABLE
      AND EXISTS "${EXP_FEMTO_3D_GSL_INCLUDE_DIR}/gsl/gsl_sf_coulomb.h"
      AND EXISTS "${EXP_FEMTO_3D_GSL_LIBRARY}"
      AND EXISTS "${EXP_FEMTO_3D_GSLCBLAS_LIBRARY}")
    set(EXP_FEMTO_3D_GSL_AVAILABLE ON)
  endif()
  if(EXISTS "${EXP_FEMTO_3D_CATS_INCLUDE_DIR}/CATS.h" AND EXISTS "${EXP_FEMTO_3D_CATS_LIBRARY}" AND EXP_FEMTO_3D_GSL_AVAILABLE)
    set(EXP_FEMTO_3D_HAS_CATS ON)
    message(STATUS "CATS finite-source Coulomb support enabled: ${EXP_FEMTO_3D_CATS_LIBRARY}")
  else()
    message(STATUS "CATS finite-source Coulomb support disabled; none/gamow modes remain available")
  endif()
endif()

add_library(exp_femto_3d_core
  src/Config.cpp
  src/InputLayout.cpp
  src/CfNormalization.cpp
  src/Cats3DModel.cpp
  src/Cats3DLambdaScan.cpp
  src/CombinedProfile.cpp
  src/Logging.cpp
  src/ProfileDisplay2D.cpp
  src/ProfileLikelihood.cpp
  src/SharedLambdaFit.cpp
  src/Workflow.cpp
)

target_include_directories(exp_femto_3d_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(exp_femto_3d_core
  PUBLIC
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::MathCore
    ROOT::Minuit
    Threads::Threads
    tomlplusplus::tomlplusplus
)

if(EXP_FEMTO_3D_HAS_CATS)
  target_compile_definitions(exp_femto_3d_core PUBLIC EXP_FEMTO_3D_HAS_CATS=1)
  target_include_directories(exp_femto_3d_core PUBLIC "${EXP_FEMTO_3D_CATS_INCLUDE_DIR}")
  if(GSL_FOUND)
    target_link_libraries(exp_femto_3d_core PUBLIC "${EXP_FEMTO_3D_CATS_LIBRARY}" GSL::gsl GSL::gslcblas)
  else()
    target_include_directories(exp_femto_3d_core PUBLIC "${EXP_FEMTO_3D_GSL_INCLUDE_DIR}")
    target_link_libraries(exp_femto_3d_core PUBLIC
      "${EXP_FEMTO_3D_CATS_LIBRARY}"
      "${EXP_FEMTO_3D_GSL_LIBRARY}"
      "${EXP_FEMTO_3D_GSLCBLAS_LIBRARY}"
    )
  endif()
  set_property(TARGET exp_femto_3d_core APPEND PROPERTY BUILD_RPATH "${EXP_FEMTO_3D_CATS_ROOT}/install/lib")
endif()

option(EXP_FEMTO_3D_ENABLE_CATS3D "Enable the local CATS3D Gaussian backend" ON)
set(EXP_FEMTO_3D_CATS3D_ROOT "${EXP_FEMTO_3D_CATS_ROOT}/CATS_3D" CACHE PATH "CATS3D helper project")
if(EXP_FEMTO_3D_ENABLE_CATS3D AND EXP_FEMTO_3D_HAS_CATS
    AND EXISTS "${EXP_FEMTO_3D_CATS_INCLUDE_DIR}/CATS3D.h"
    AND EXISTS "${EXP_FEMTO_3D_CATS3D_ROOT}/CATS3DFastKernelCache.cpp")
  find_package(ROOT QUIET COMPONENTS Physics)
  if(TARGET ROOT::Physics)
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_INCLUDES "${EXP_FEMTO_3D_CATS_INCLUDE_DIR}")
  set(CMAKE_REQUIRED_LIBRARIES "${EXP_FEMTO_3D_CATS_LIBRARY}" ROOT::Physics ROOT::Core)
  unset(EXP_FEMTO_3D_CATS3D_API_MATCHES CACHE)
  check_cxx_source_compiles("#include <CATS3D.h>
    int main() { CATS3D c; c.SetFastKernelSubdivisions(1); c.SetAnaSource([](void*, double*) { return 1.0; }, &c, 6); c.BuildFastKernelOnly(); return c.GetFastWaveFunctionKernel().empty(); }"
    EXP_FEMTO_3D_CATS3D_API_MATCHES)
  unset(CMAKE_REQUIRED_INCLUDES)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(EXP_FEMTO_3D_CATS3D_API_MATCHES)
    target_sources(exp_femto_3d_core PRIVATE "${EXP_FEMTO_3D_CATS3D_ROOT}/CATS3DFastKernelCache.cpp")
    target_include_directories(exp_femto_3d_core PRIVATE "${EXP_FEMTO_3D_CATS3D_ROOT}")
    target_link_libraries(exp_femto_3d_core PUBLIC ROOT::Physics)
    target_compile_definitions(exp_femto_3d_core PUBLIC EXP_FEMTO_3D_HAS_CATS3D=1)
    target_compile_definitions(exp_femto_3d_core PRIVATE
      EXP_FEMTO_3D_CATS3D_CACHE_DIR="${CMAKE_CURRENT_BINARY_DIR}/cats3d_kernels")
    message(STATUS "CATS3D Gaussian backend enabled")
  endif()
  endif()
endif()

target_compile_options(exp_femto_3d_core
  PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

target_compile_definitions(exp_femto_3d_core
  PRIVATE EXP_FEMTO_3D_EXECUTABLE_PATH="${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/exp_femto_3d"
)

add_executable(exp_femto_3d
  app/main.cpp
)

target_link_libraries(exp_femto_3d
  PRIVATE
    exp_femto_3d_core
)

# Offline rendering and harmonic summaries; never reruns correlation-function fits.
add_executable(cats3d_scan_report app/Cats3DScanReport.cpp)
target_include_directories(cats3d_scan_report PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(cats3d_scan_report PRIVATE exp_femto_3d_core)


add_executable(coulomb_kernel_validation_test
  tests/CoulombKernelValidationTest.cpp
)

target_link_libraries(coulomb_kernel_validation_test
  PRIVATE
    exp_femto_3d_core
)

add_test(
  NAME coulomb_kernel_validation_test
  COMMAND coulomb_kernel_validation_test
)

add_executable(config_parse_validation_test
  tests/ConfigParseValidationTest.cpp
)

target_link_libraries(config_parse_validation_test
  PRIVATE
    exp_femto_3d_core
)

target_compile_definitions(config_parse_validation_test
  PRIVATE EXP_FEMTO_3D_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

add_test(
  NAME config_parse_validation_test
  COMMAND config_parse_validation_test
)

add_executable(profile_likelihood_driver_test
  tests/ProfileLikelihoodDriverTest.cpp
  src/ProfileLikelihood.cpp
)

target_include_directories(profile_likelihood_driver_test
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_test(
  NAME profile_likelihood_driver_test
  COMMAND profile_likelihood_driver_test
)

add_executable(profile_display_2d_test
  tests/ProfileDisplay2DTest.cpp
  src/ProfileDisplay2D.cpp
)

target_include_directories(profile_display_2d_test
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_test(
  NAME profile_display_2d_test
  COMMAND profile_display_2d_test
)

add_executable(shared_lambda_fit_test
  tests/SharedLambdaFitTest.cpp
  src/SharedLambdaFit.cpp
)

target_include_directories(shared_lambda_fit_test
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_test(
  NAME shared_lambda_fit_test
  COMMAND shared_lambda_fit_test
)

add_executable(combined_profile_driver_test
  tests/CombinedProfileDriverTest.cpp
  src/CombinedProfile.cpp
  src/ProfileLikelihood.cpp
)

target_include_directories(combined_profile_driver_test
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_test(
  NAME combined_profile_driver_test
  COMMAND combined_profile_driver_test
)

add_executable(progress_render_test
  tests/ProgressRenderTest.cpp
)

target_link_libraries(progress_render_test
  PRIVATE
    exp_femto_3d_core
)

add_test(
  NAME progress_render_test
  COMMAND progress_render_test
)

add_executable(root_runtime_probe
  tests/RootRuntimeProbe.cpp
)

target_link_libraries(root_runtime_probe
  PRIVATE
    exp_femto_3d_core
)

add_executable(slice_catalog_roundtrip_test
  tests/SliceCatalogRoundTripTest.cpp
)

target_link_libraries(slice_catalog_roundtrip_test
  PRIVATE
    exp_femto_3d_core
)

add_test(
  NAME slice_catalog_roundtrip_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe>
          $<TARGET_FILE:slice_catalog_roundtrip_test>
)

add_executable(build_cf_rebin_test
  tests/BuildCfRebinTest.cpp
)

target_link_libraries(build_cf_rebin_test
  PRIVATE
    exp_femto_3d_core
)

add_test(
  NAME build_cf_rebin_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe>
          $<TARGET_FILE:build_cf_rebin_test>
)

add_executable(input_format_test tests/InputFormatTest.cpp)
target_link_libraries(input_format_test PRIVATE exp_femto_3d_core)
add_test(NAME input_format_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe> $<TARGET_FILE:input_format_test>)
set_tests_properties(input_format_test PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 180)

add_executable(workflow_smoke_test
  tests/WorkflowSmokeTest.cpp
)

target_link_libraries(workflow_smoke_test
  PRIVATE
    exp_femto_3d_core
)

add_test(
  NAME workflow_smoke_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe>
          $<TARGET_FILE:workflow_smoke_test>
)

set_tests_properties(
  build_cf_rebin_test
  slice_catalog_roundtrip_test
  workflow_smoke_test
  PROPERTIES
    SKIP_RETURN_CODE 77
)

add_executable(cats3d_scan_test tests/Cats3DLambdaScanTest.cpp src/Cats3DLambdaScan.cpp src/ProfileLikelihood.cpp)
target_include_directories(cats3d_scan_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
add_test(NAME cats3d_scan_test COMMAND cats3d_scan_test)
add_executable(cats3d_model_test tests/Cats3DModelTest.cpp)
target_include_directories(cats3d_model_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(cats3d_model_test PRIVATE exp_femto_3d_core)
add_test(NAME cats3d_model_test COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh"
  $<TARGET_FILE:root_runtime_probe> $<TARGET_FILE:cats3d_model_test>)
set_tests_properties(cats3d_model_test PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 180)

add_executable(cats3d_workflow_test tests/Cats3DWorkflowTest.cpp)
target_include_directories(cats3d_workflow_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_compile_definitions(cats3d_workflow_test PRIVATE EXP_FEMTO_3D_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(cats3d_workflow_test PRIVATE exp_femto_3d_core)
add_test(NAME cats3d_workflow_test COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh"
  $<TARGET_FILE:root_runtime_probe> $<TARGET_FILE:cats3d_workflow_test>)
set_tests_properties(cats3d_workflow_test PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600)

# Shell normalization depends on both ROOT histogram arithmetic and the count-likelihood workflow.
add_executable(cf_normalization_test tests/CfNormalizationTest.cpp)
target_link_libraries(cf_normalization_test PRIVATE exp_femto_3d_core)
add_test(NAME cf_normalization_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe> $<TARGET_FILE:cf_normalization_test>)
set_tests_properties(cf_normalization_test PROPERTIES SKIP_RETURN_CODE 77)

cmake_minimum_required(VERSION 3.16)

project(ExpFemto3D LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json" FORCE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

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
  src/Logging.cpp
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

target_compile_options(exp_femto_3d_core
  PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

add_executable(exp_femto_3d
  app/main.cpp
)

target_link_libraries(exp_femto_3d
  PRIVATE
    exp_femto_3d_core
)


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

add_test(
  NAME config_parse_validation_test
  COMMAND config_parse_validation_test
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

cmake_minimum_required(VERSION 3.16)

project(ExpFemto1D LANGUAGES CXX)
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

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad Tree TreePlayer MathCore Minuit)
find_package(GSL REQUIRED)
find_package(tomlplusplus CONFIG REQUIRED
  HINTS
    /opt/homebrew/lib/cmake/tomlplusplus
    /opt/homebrew/opt/tomlplusplus
)

set(EXP_FEMTO_1D_DEFAULT_CATS_ROOT "/Users/allenzhou/Research_software/CATS/install")
set(CATS_ROOT "${CATS_ROOT}" CACHE PATH "Installed CATS prefix")
if(NOT CATS_ROOT)
  set(CATS_ROOT "${EXP_FEMTO_1D_DEFAULT_CATS_ROOT}")
endif()

find_path(CATS_INCLUDE_DIR
  NAMES CATS.h
  HINTS
    "${CATS_ROOT}/include"
    "${EXP_FEMTO_1D_DEFAULT_CATS_ROOT}/include"
)

find_library(CATS_LIBRARY
  NAMES CATS libCATS.dylib
  HINTS
    "${CATS_ROOT}/lib"
    "${EXP_FEMTO_1D_DEFAULT_CATS_ROOT}/lib"
)

if(NOT CATS_INCLUDE_DIR OR NOT CATS_LIBRARY)
  message(FATAL_ERROR "Failed to locate installed CATS headers or library. Set CATS_ROOT to the install prefix.")
endif()

get_filename_component(CATS_LIBRARY_DIR "${CATS_LIBRARY}" DIRECTORY)
list(APPEND CMAKE_BUILD_RPATH "${CATS_LIBRARY_DIR}")

add_library(exp_femto_1d_core
  src/CatsModel.cpp
  src/Config.cpp
  src/Logging.cpp
  src/Workflow.cpp
)

target_include_directories(exp_femto_1d_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
    ${CATS_INCLUDE_DIR}
)

target_link_libraries(exp_femto_1d_core
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
    tomlplusplus::tomlplusplus
    GSL::gsl
    GSL::gslcblas
    ${CATS_LIBRARY}
)

target_compile_options(exp_femto_1d_core
  PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

add_executable(exp_femto_1d
  app/main.cpp
)

target_link_libraries(exp_femto_1d
  PRIVATE
    exp_femto_1d_core
)

add_executable(config_parse_validation_test
  tests/ConfigParseValidationTest.cpp
)

target_link_libraries(config_parse_validation_test
  PRIVATE
    exp_femto_1d_core
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
    exp_femto_1d_core
)

add_executable(slice_catalog_roundtrip_test
  tests/SliceCatalogRoundTripTest.cpp
)

target_link_libraries(slice_catalog_roundtrip_test
  PRIVATE
    exp_femto_1d_core
)

add_test(
  NAME slice_catalog_roundtrip_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe>
          $<TARGET_FILE:slice_catalog_roundtrip_test>
)

add_executable(workflow_smoke_test
  tests/WorkflowSmokeTest.cpp
)

target_link_libraries(workflow_smoke_test
  PRIVATE
    exp_femto_1d_core
)

add_test(
  NAME workflow_smoke_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe>
          $<TARGET_FILE:workflow_smoke_test>
)

add_executable(cats_fit_smoke_test
  tests/CatsFitSmokeTest.cpp
)

target_link_libraries(cats_fit_smoke_test
  PRIVATE
    exp_femto_1d_core
)

add_test(
  NAME cats_fit_smoke_test
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_root_guard.sh
          $<TARGET_FILE:root_runtime_probe>
          $<TARGET_FILE:cats_fit_smoke_test>
)

set_tests_properties(
  slice_catalog_roundtrip_test
  workflow_smoke_test
  cats_fit_smoke_test
  PROPERTIES
    SKIP_RETURN_CODE 77
)

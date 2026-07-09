cmake_minimum_required(VERSION 3.20)

project(geodesic_draping
  VERSION 0.1.0
  LANGUAGES CXX)

option(GEODESIC_DRAPING_BUILD_TESTS "Build geodesic draping tests" ON)
option(GEODESIC_DRAPING_ENABLE_POLYSCOPE "Enable optional Polyscope plotting helpers" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(extern/geometry-central)

if(GEODESIC_DRAPING_ENABLE_POLYSCOPE)
  set(POLYSCOPE_BACKEND_OPENGL3_EGL OFF CACHE STRING "Disable EGL backend on Windows" FORCE)
  add_subdirectory(extern/polyscope)
endif()

add_library(geodesic_draping
  src/diagnostics.cpp
  src/field_processing.cpp
  src/geodrape.cpp
  src/geometrycentral_adapter.cpp
  src/generator_tracing.cpp
  src/plotting.cpp
  src/quality.cpp
  src/seed_projection.cpp
  src/signed_heat.cpp
  src/custom_signed_heat.cpp
  src/version.cpp)

target_include_directories(geodesic_draping
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(geodesic_draping
  PUBLIC
    geometry-central)

if(GEODESIC_DRAPING_ENABLE_POLYSCOPE)
  target_link_libraries(geodesic_draping
    PUBLIC
      polyscope)
  target_compile_definitions(geodesic_draping
    PUBLIC
      GEODESIC_DRAPING_HAS_POLYSCOPE=1)
else()
  target_compile_definitions(geodesic_draping
    PUBLIC
      GEODESIC_DRAPING_HAS_POLYSCOPE=0)
endif()

if(GEODESIC_DRAPING_BUILD_TESTS)
  enable_testing()

  add_executable(smoke_eigen
    tests/smoke_eigen.cpp)
  target_link_libraries(smoke_eigen
    PRIVATE
      geodesic_draping)
  add_test(NAME smoke_eigen COMMAND smoke_eigen)

  add_executable(smoke_geometrycentral
    tests/smoke_geometrycentral.cpp)
  target_link_libraries(smoke_geometrycentral
    PRIVATE
      geodesic_draping)
  add_test(NAME smoke_geometrycentral COMMAND smoke_geometrycentral)

  add_executable(test_seed_projection
    tests/test_seed_projection.cpp)
  target_link_libraries(test_seed_projection
    PRIVATE
      geodesic_draping)
  target_include_directories(test_seed_projection
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(test_seed_projection
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
  add_test(NAME test_seed_projection COMMAND test_seed_projection)

  add_executable(test_geometrycentral_adapter
    tests/test_geometrycentral_adapter.cpp)
  target_link_libraries(test_geometrycentral_adapter
    PRIVATE
      geodesic_draping)
  target_include_directories(test_geometrycentral_adapter
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(test_geometrycentral_adapter
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
  add_test(NAME test_geometrycentral_adapter COMMAND test_geometrycentral_adapter)

  add_executable(test_generator_tracing
    tests/test_generator_tracing.cpp)
  target_link_libraries(test_generator_tracing
    PRIVATE
      geodesic_draping)
  target_include_directories(test_generator_tracing
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(test_generator_tracing
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
  add_test(NAME test_generator_tracing COMMAND test_generator_tracing)

  add_executable(test_signed_heat
    tests/test_signed_heat.cpp)
  target_link_libraries(test_signed_heat
    PRIVATE
      geodesic_draping)
  target_include_directories(test_signed_heat
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(test_signed_heat
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
  add_test(NAME test_signed_heat COMMAND test_signed_heat)

  add_executable(test_custom_signed_heat
    tests/test_custom_signed_heat.cpp)
  target_link_libraries(test_custom_signed_heat
    PRIVATE
      geodesic_draping)
  target_include_directories(test_custom_signed_heat
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(test_custom_signed_heat
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
  add_test(NAME test_custom_signed_heat COMMAND test_custom_signed_heat)

  add_executable(test_field_processing
    tests/test_field_processing.cpp)
  target_link_libraries(test_field_processing
    PRIVATE
      geodesic_draping)
  target_include_directories(test_field_processing
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(test_field_processing
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
  add_test(NAME test_field_processing COMMAND test_field_processing)

  add_executable(test_diagnostics
    tests/test_diagnostics.cpp)
  target_link_libraries(test_diagnostics
    PRIVATE
      geodesic_draping)
  target_include_directories(test_diagnostics
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(test_diagnostics
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
  add_test(NAME test_diagnostics COMMAND test_diagnostics)

  add_executable(test_quality
    tests/test_quality.cpp)
  target_link_libraries(test_quality
    PRIVATE
      geodesic_draping)
  target_include_directories(test_quality
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(test_quality
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
  add_test(NAME test_quality COMMAND test_quality)

  add_executable(test_complete_drape
    tests/test_complete_drape.cpp)
  target_link_libraries(test_complete_drape
    PRIVATE
      geodesic_draping)
  target_include_directories(test_complete_drape
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(test_complete_drape
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
  add_test(NAME test_complete_drape COMMAND test_complete_drape)
endif()

if(GEODESIC_DRAPING_ENABLE_POLYSCOPE)
  add_executable(debug_geodrape
    tools/debug_geodrape.cpp)
  target_link_libraries(debug_geodrape
    PRIVATE
      geodesic_draping)
  target_include_directories(debug_geodrape
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  target_compile_definitions(debug_geodrape
    PRIVATE
      GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")
endif()

add_executable(benchmark_geodrape
  tools/benchmark_geodrape.cpp)
target_link_libraries(benchmark_geodrape
  PRIVATE
    geodesic_draping)
target_include_directories(benchmark_geodrape
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tools)
target_compile_definitions(benchmark_geodrape
  PRIVATE
    GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")

add_executable(profile_geodrape
  tools/profile_geodrape.cpp)
target_link_libraries(profile_geodrape
  PRIVATE
    geodesic_draping)
target_include_directories(profile_geodrape
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tools)
target_compile_definitions(profile_geodrape
  PRIVATE
    GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")

add_executable(debug_intrinsic_subdivision
  tools/debug_intrinsic_subdivision.cpp)
target_link_libraries(debug_intrinsic_subdivision
  PRIVATE
    geodesic_draping)
target_include_directories(debug_intrinsic_subdivision
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tools)
target_compile_definitions(debug_intrinsic_subdivision
  PRIVATE
    GEODESIC_DRAPING_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/fixtures")

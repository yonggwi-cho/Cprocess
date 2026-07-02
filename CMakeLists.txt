cmake_minimum_required(VERSION 3.16)
project(cprocess LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP)
find_package(Threads REQUIRED)

add_library(cprocess_core
  src/mesh.cpp
  src/box_mesh.cpp
  src/cell_locator.cpp
  src/gmsh_reader.cpp
  src/sparse.cpp
  src/topology.cpp
  src/remesh.cpp
  src/ale_mover.cpp
  src/oxidation.cpp
  src/field_transfer.cpp
  src/levelset.cpp
  src/mechanics.cpp
  src/materials.cpp
  src/implant.cpp
  src/mc_implant.cpp
  src/diffusion.cpp
  src/vtk_writer.cpp
  src/process.cpp
  src/deck.cpp)
target_include_directories(cprocess_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(cprocess_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(cprocess_core PUBLIC Threads::Threads)
if(OpenMP_CXX_FOUND)
  target_link_libraries(cprocess_core PUBLIC OpenMP::OpenMP_CXX)
  message(STATUS "OpenMP enabled")
else()
  message(STATUS "OpenMP not found — running single-threaded")
endif()
if(NOT MSVC)
  target_compile_options(cprocess_core PRIVATE -Wall -Wextra)
endif()

add_executable(cprocess src/main.cpp)
target_link_libraries(cprocess PRIVATE cprocess_core)

add_executable(nmos_locos_flow examples/nmos_locos_flow.cpp)
target_link_libraries(nmos_locos_flow PRIVATE cprocess_core)

enable_testing()
foreach(t mesh solver gmsh diffusion mc remesh ale oxidation field_transfer levelset mechanics resist photo integration ted units flow oxidize_flow segregation multimaterial activation rta pearson)
  add_executable(test_${t} tests/test_${t}.cpp)
  target_link_libraries(test_${t} PRIVATE cprocess_core)
  add_test(NAME ${t} COMMAND test_${t})
endforeach()

# Python bindings (optional — requires pybind11)
find_package(pybind11 CONFIG QUIET)
if(pybind11_FOUND)
  pybind11_add_module(_cprocess python/_cprocess.cpp)
  target_link_libraries(_cprocess PRIVATE cprocess_core)
  # Install the .so next to the cprocess Python package directory
  install(TARGETS _cprocess
          LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/python/cprocess)
  message(STATUS "Python bindings enabled")
else()
  message(STATUS "pybind11 not found — Python bindings skipped")
endif()

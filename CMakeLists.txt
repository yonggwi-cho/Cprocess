cmake_minimum_required(VERSION 3.16)
project(cprocess LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP)
find_package(Threads REQUIRED)

option(CPROCESS_GPU "Enable OpenMP target offload" OFF)
option(CPROCESS_MPI "Enable MPI-distributed MC implant" OFF)

add_library(cprocess_core
  src/mesh.cpp
  src/box_mesh.cpp
  src/cell_locator.cpp
  src/gmsh_reader.cpp
  src/gds_reader.cpp
  src/sparse.cpp
  src/topology.cpp
  src/remesh.cpp
  src/ale_mover.cpp
  src/oxidation.cpp
  src/field_transfer.cpp
  src/levelset.cpp
  src/mechanics.cpp
  src/fem.cpp
  src/materials.cpp
  src/param_db.cpp
  src/implant.cpp
  src/mc_implant.cpp
  src/diffusion.cpp
  src/oxidant_solver.cpp
  src/vtk_writer.cpp
  src/state_io.cpp
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
if(CPROCESS_GPU)
  target_compile_definitions(cprocess_core PUBLIC CPROCESS_GPU)
  # gcc: -foffload=nvptx-none / clang: -fopenmp-targets=nvptx64
  # Detection is not attempted here; CPROCESS_OFFLOAD_FLAGS is a
  # user-supplied cache variable naming the offload flags for the
  # local toolchain.
  set(CPROCESS_OFFLOAD_FLAGS "" CACHE STRING "OpenMP target offload flags")
  target_compile_options(cprocess_core PUBLIC ${CPROCESS_OFFLOAD_FLAGS})
  target_link_options(cprocess_core PUBLIC ${CPROCESS_OFFLOAD_FLAGS})
endif()
if(CPROCESS_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
  target_link_libraries(cprocess_core PUBLIC MPI::MPI_CXX)
  target_compile_definitions(cprocess_core PUBLIC CPROCESS_MPI)
endif()

add_executable(cprocess src/main.cpp)
target_link_libraries(cprocess PRIVATE cprocess_core)

add_executable(nmos_locos_flow examples/nmos_locos_flow.cpp)
target_link_libraries(nmos_locos_flow PRIVATE cprocess_core)

enable_testing()
foreach(t mesh solver gmsh diffusion mc compound_bca remesh ale oxidation field_transfer levelset mechanics resist photo integration ted units flow oxidize_flow segregation multimaterial activation rta pearson params etch_depo state_io gds dopant_clusters new_dopants oed oxidant ox2d fem topo amg block_solver topology_parallel newton adaptive_dt species_parallel rcm)
  add_executable(test_${t} tests/test_${t}.cpp)
  target_link_libraries(test_${t} PRIVATE cprocess_core)
  add_test(NAME ${t} COMMAND test_${t})
endforeach()

if(CPROCESS_GPU)
  foreach(t gpu)
    add_executable(test_${t} tests/test_${t}.cpp)
    target_link_libraries(test_${t} PRIVATE cprocess_core)
    add_test(NAME ${t} COMMAND test_${t})
  endforeach()
endif()

if(CPROCESS_MPI)
  add_executable(test_mpi_mc tests/test_mpi_mc.cpp)
  target_link_libraries(test_mpi_mc PRIVATE cprocess_core)
  add_test(NAME mpi_mc COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
           ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_mc> ${MPIEXEC_POSTFLAGS})
endif()

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

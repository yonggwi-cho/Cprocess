cmake_minimum_required(VERSION 3.16)
project(cprocess LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(cprocess_core
  src/mesh.cpp
  src/box_mesh.cpp
  src/gmsh_reader.cpp
  src/sparse.cpp
  src/materials.cpp
  src/implant.cpp
  src/diffusion.cpp
  src/vtk_writer.cpp
  src/deck.cpp)
target_include_directories(cprocess_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(NOT MSVC)
  target_compile_options(cprocess_core PRIVATE -Wall -Wextra)
endif()

add_executable(cprocess src/main.cpp)
target_link_libraries(cprocess PRIVATE cprocess_core)

enable_testing()
foreach(t mesh solver gmsh diffusion)
  add_executable(test_${t} tests/test_${t}.cpp)
  target_link_libraries(test_${t} PRIVATE cprocess_core)
  add_test(NAME ${t} COMMAND test_${t})
endforeach()

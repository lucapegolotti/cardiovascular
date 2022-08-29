
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 11)

PROJECT(slice-results)

find_package(VTK REQUIRED PATHS /home/lpegolotti/VTK-install/lib64/cmake/vtk-9.1/)
include(${VTK_USE_FILE})

add_executable(slice-results
    main_slice_results.cpp
    Centerlines.h Centerlines.cpp
    Graphics.h Graphics.cpp
    Mesh.h Mesh.cpp
    Slice.h Slice.cpp
)

add_executable(extract-all-slices
    main_extract_all_slices.cpp
    Centerlines.h Centerlines.cpp
    Graphics.h Graphics.cpp
    Mesh.h Mesh.cpp
    Slice.h Slice.cpp
)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp")
target_link_libraries(slice-results ${VTK_LIBRARIES})
target_link_libraries(extract-all-slices ${VTK_LIBRARIES})

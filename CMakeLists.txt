cmake_minimum_required(VERSION 3.16)
project(3D-Placement)

set(CMAKE_CXX_STANDARD 17)
set(SOURCES
        src/main.cpp

        src/Chip.cpp
        src/dataStructures.cpp
        src/Partitioner.cpp

        src/tool/fft.cpp
        src/tool/fftsg.cpp
        src/tool/fftsg2d.cpp
        src/tool/Drawer.cpp
        )

set(INCLUDE_DIRS
        ${PROJECT_SOURCE_DIR}/include/algorithm
        ${PROJECT_SOURCE_DIR}/include/circuit
        ${PROJECT_SOURCE_DIR}/include/dataStructures
        ${PROJECT_SOURCE_DIR}/include/tool
        ${PROJECT_SOURCE_DIR}/submodules/CImg
        )

add_executable(placer ${SOURCES})
target_include_directories(placer PRIVATE ${INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
add_subdirectory(submodules/OpenROAD)
target_link_libraries(placer
        PUBLIC
        odb
        mpl2_lib
        Eigen3::Eigen
        X11
        )
target_include_directories(mpl2_lib
        PUBLIC
        submodules/OpenROAD/src/mpl2/src
        )
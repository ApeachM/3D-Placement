cmake_minimum_required(VERSION 3.16)
project(3D-Placement)

set(CMAKE_CXX_STANDARD 17)
set(SOURCES
        src/main.cpp

        src/circuit/Chip.cpp
        src/circuit/Parser.cpp

        src/dataStructures/Die.cpp
        src/dataStructures/Instance.cpp
        src/dataStructures/Net.cpp
        src/dataStructures/Pin.cpp
        src/dataStructures/Bin.cpp

        src/algorithm/do3DPlace.cpp
        src/algorithm/normalPlace.cpp
        src/algorithm/InitialPlacer.cpp
        src/algorithm/NesterovPlacer.cpp

        src/tool/fft.cpp
        src/tool/fftsg.cpp
        src/tool/fftsg2d.cpp

        src/etc/dbTutorial.cpp
        )

set(INCLUDE_DIRS
        ${PROJECT_SOURCE_DIR}/include/algorithm
        ${PROJECT_SOURCE_DIR}/include/circuit
        ${PROJECT_SOURCE_DIR}/include/dataStructures
        ${PROJECT_SOURCE_DIR}/include/tool
        )

add_executable(placer ${SOURCES})
target_include_directories(placer PRIVATE ${INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
add_subdirectory(submodule/OpenROAD)
target_link_libraries(placer
        PUBLIC
#        openroad
        odb
#        mpl2
        Eigen3::Eigen
        )

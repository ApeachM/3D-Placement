cmake_minimum_required(VERSION 3.14..3.16)
project(3D-Placement)

set(CMAKE_CXX_STANDARD 11)
set(SOURCES
        src/main.cpp

        src/circuit/Chip.cpp
        src/circuit/Parser.cpp

        src/dataStructures/Die.cpp
        src/dataStructures/Instance.cpp
        src/dataStructures/Net.cpp
        src/dataStructures/Pin.cpp

        src/algorithm/place.cpp

        src/etc/dbTutorial.cpp
        src/algorithm/normalPlace.cpp)

set(INCLUDE_DIRS
        #${PROJECT_SOURCE_DIR}/include
        #${PROJECT_SOURCE_DIR}/include/algorithm
        ${PROJECT_SOURCE_DIR}/include/circuit
        ${PROJECT_SOURCE_DIR}/include/dataStructures
        ${PROJECT_SOURCE_DIR}/submodule/CImg
        )

add_executable(placer ${SOURCES})
target_include_directories(placer PRIVATE ${INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
add_subdirectory(submodule/OpenDB)
add_subdirectory(submodule/igraph)
target_link_libraries(placer
        PRIVATE
        Eigen3::Eigen
        opendb
        X11
        pthread
        igraph
        )

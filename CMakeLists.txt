cmake_minimum_required(VERSION 3.14..3.16)
project(3D-Placement-with-D2D-Vertical-Connections)

set(CMAKE_CXX_STANDARD 11)
set(SOURCES
        src/main.cpp

        src/circuit/Circuit.cpp
        src/circuit/Parser.cpp

        src/dataStructures/Die.cpp
        src/dataStructures/Instance.cpp
        src/dataStructures/Net.cpp
        src/dataStructures/Pin.cpp

        src/algorithm/place.cpp

        src/etc/dbTutorial.cpp
        src/dataStructures/D2DChip.cpp
        )

set(INCLUDE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/include/algorithm
        ${PROJECT_SOURCE_DIR}/include/circuit
        ${PROJECT_SOURCE_DIR}/include/dataStructures
        ${PROJECT_SOURCE_DIR}/submodule/CImg
        )

add_executable(placer ${SOURCES})
target_include_directories(placer PRIVATE ${INCLUDE_DIRS})

add_subdirectory(submodule/OpenDB)
target_link_libraries(placer PRIVATE opendb X11 pthread)

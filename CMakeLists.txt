cmake_minimum_required(VERSION 3.14..3.16)
project(placer)

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
        )

add_executable(placer ${SOURCES})
target_include_directories(placer PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(placer PRIVATE ${PROJECT_SOURCE_DIR}/include/circuit)
target_include_directories(placer PRIVATE ${PROJECT_SOURCE_DIR}/include/dataStructures)

add_subdirectory(submodule/OpenDB)
target_link_libraries(placer PRIVATE opendb)

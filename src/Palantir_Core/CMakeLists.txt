cmake_minimum_required (VERSION 3.0)

project(Palantir_Core)

set(CMAKE_CXX_STANDARD 11)
find_package(Armadillo REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})

file(GLOB HEADERS *.hpp)
file(GLOB SOURCES *.cpp)
file(GLOB TEST_SOURCES test/*.cpp test/*.hpp)

add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE "include/")

add_executable(test-Palantir ${HEADERS} ${TEST_SOURCES})
add_library(Palantir SHARED ${SOURCES} ${HEADERS})

target_link_libraries(Palantir ${ARMADILLO_LIBRARIES})
target_link_libraries(test-Palantir ${ARMADILLO_LIBRARIES} Palantir Catch)

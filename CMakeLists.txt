cmake_minimum_required(VERSION 3.0)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/"
    CACHE INTERNAL "Location of our custom CMake modules.")
set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(AOS)

add_subdirectory(code)
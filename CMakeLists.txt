cmake_minimum_required(VERSION 2.8)
project(PPData CXX)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default Build Type is Release" FORCE)

# build unittest
include_directories(src/
                    3rdparty/googletest-release-1.7.0/include/)
set(gtest_disable_pthreads ON CACHE BOOL "disable pthreads")
add_subdirectory(3rdparty/googletest-release-1.7.0)

add_definitions("-std=c++1y")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

if (MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
endif()

add_executable(unittest test/Test_PPData.cpp src/PPData.cpp)
target_link_libraries(unittest gtest gtest_main)

set(CMAKE_DEBUG_POSTFIX "d")
add_library(ppdata STATIC src/PPData.cpp)
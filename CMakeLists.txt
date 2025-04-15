cmake_minimum_required(VERSION 3.16)
project(BFL)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-pipe -march=native -Wfatal-errors -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

add_compile_options(-Wno-gnu-anonymous-struct -Wno-nested-anon-types)

set(CMAKE_CXX_FLAGS_RELEASE "-ffast-math -DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-Werror -Wpedantic -pedantic-errors -Wall -Wextra -g -fno-omit-frame-pointer -O0")
set(CMAKE_CXX_FLAGS_SANITIZE "-fsanitize=address,undefined -fno-omit-frame-pointer -g -O1")
set(CMAKE_LINKER_FLAGS_SANITIZE "-fsanitize=address,undefined")

if(CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_SANITIZE}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS_SANITIZE}")
    message(STATUS "Sanitizer build enabled")
endif()

add_executable(BFL main.cpp)
target_include_directories(BFL PUBLIC datastructures external)
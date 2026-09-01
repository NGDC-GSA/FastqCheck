# CMakeLists.txt of FastqCheck
#
# usage:
#   cmake -S . -B build && cmake --build build -j
#       configure + build fastqcheck (default: Release)

cmake_minimum_required(VERSION 3.16)
project(FastqCheck C)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_C_STANDARD 99)          # same as the old Makefile's -std=c99
set(CMAKE_C_STANDARD_REQUIRED ON)

# ---- zlib-ng: build the local source tree as a static library -------------
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Build test binaries" FORCE)
set(ZLIB_COMPAT ON CACHE BOOL "Compile with zlib compatible API" FORCE)
set(WITH_GZFILEOP ON CACHE BOOL "Compile with support for gzFile related functions" FORCE)
add_subdirectory(external/zlibng) # target: zlib-ng::zlib

# ---- bzip2: build the local source tree as a static library ---------------
set(ENABLE_LIB_ONLY ON CACHE BOOL "Build libbz2 only" FORCE)
set(ENABLE_STATIC_LIB ON CACHE BOOL "Build libbz2 in static mode also" FORCE)
set(ENABLE_SHARED_LIB OFF CACHE BOOL "Build libbz2 as a shared library" FORCE)
set(ENABLE_TESTS OFF CACHE BOOL "Build/enable unit tests" FORCE)
set(ENABLE_DOCS OFF CACHE BOOL "Generate documentation" FORCE)
add_subdirectory(external/bzip2)  # target: bz2_static

# disable the warnings of bzip2 during compile
if(MSVC)
    target_compile_options(bz2_static PRIVATE /w)
else()
    target_compile_options(bz2_static PRIVATE -w)
endif()


# ---- system libraries ------------------------------------------------------
find_package(Threads REQUIRED)

# ---- fastqcheck -----------------------------------------------------------
add_executable(fastqcheck
    bloom_filter.c
    file_read.c
    kqueue.c
    params_parse.c
    threadpool.c
    fastq_check.c)

target_link_libraries(fastqcheck PRIVATE
    zlib-ng::zlib
    bz2_static
    m
    Threads::Threads)
# bz2_static does not export include dirs, so point <bzlib.h> at the
# vendored copy's public header.
target_include_directories(fastqcheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/bzip2)
target_compile_definitions(fastqcheck PRIVATE _GNU_SOURCE)

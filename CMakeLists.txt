# CMakeLists.txt of FastqCheck
#
# usage:
#   cmake -S . -B build && cmake --build build -j
#       configure + build fastq_check (default: Release)

cmake_minimum_required(VERSION 3.16)
project(FastqCheck C)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_C_STANDARD 99)          # same as the old Makefile's -std=c99
set(CMAKE_C_STANDARD_REQUIRED ON)

# ---- zlib-ng: build the local source tree as a static library -------------
# ZLIB_COMPAT=ON makes zlib-ng export the classic zlib API (gzopen/gzread/
# gzwrite/gzclose/... under their original names, with a generated zlib.h),
# so FastqCheck's existing zlib calls compile and link without source changes.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Build test binaries" FORCE)
set(ZLIB_COMPAT ON CACHE BOOL "Compile with zlib compatible API" FORCE)
set(WITH_GZFILEOP ON CACHE BOOL "Compile with support for gzFile related functions" FORCE)
add_subdirectory(zlibng)          # target: zlib-ng::zlib

# ---- system libraries ------------------------------------------------------
find_package(BZip2 REQUIRED)
find_package(Threads REQUIRED)

# ---- fastq_check -----------------------------------------------------------
add_executable(fastq_check
    bloom_filter.c
    file_read.c
    kqueue.c
    params_parse.c
    threadpool.c
    fastq_check.c)

target_link_libraries(fastq_check PRIVATE
    zlib-ng::zlib
    BZip2::BZip2
    m
    Threads::Threads)
target_compile_definitions(fastq_check PRIVATE _GNU_SOURCE)

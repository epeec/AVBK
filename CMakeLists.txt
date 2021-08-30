cmake_minimum_required(VERSION 3.5)
project (gs_ompss)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set(HDF5_USE_STATIC_LIBRARIES TRUE)
find_package (HDF5 1.8.11 COMPONENTS C REQUIRED)
find_package (METIS REQUIRED)
find_package(OpenMP)

include_directories (.)

## ASAN settings
include (CheckCXXCompilerFlag)
set (CMAKE_REQUIRED_FLAGS "-fsanitize=address -fsanitize=undefined -fsanitize=leak -fno-omit-frame-pointer")
CHECK_CXX_COMPILER_FLAG ("-fsanitize=address -fsanitize=undefined -fsanitize=leak -fno-omit-frame-pointer" HAVE_FLAG_SANITIZE_ADDRESS)
unset (CMAKE_REQUIRED_FLAGS)
if (HAVE_FLAG_SANITIZE_ADDRESS)
  set (ASAN_FLAGS "-fsanitize=address -fsanitize=undefined -fsanitize=leak -fno-omit-frame-pointer")
else ()
  set(ASAN_FLAGS)
  message(WARNING "ASAN not found cannot check for memory leaks in debug mode!")
endif (HAVE_FLAG_SANITIZE_ADDRESS)

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ASAN_FLAGS} ${OpenMP_C_FLAGS} -O2")

add_executable(gs_ompss
   backend.cpp
   kernel.h
   main.cpp)

target_include_directories (gs_ompss PUBLIC
   ${HDF5_INCLUDE_DIR}
   ${METIS_INCLUDE_DIR})

target_link_libraries (gs_ompss PUBLIC
   ${HDF5_LIBRARIES}
   ${METIS_LIBRARY}
   ${CMAKE_DL_LIBS})

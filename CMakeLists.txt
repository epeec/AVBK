cmake_minimum_required(VERSION 3.5)
project (AVBP-mini-kernel)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set(HDF5_USE_STATIC_LIBRARIES TRUE)
find_package (HDF5 1.8.11 COMPONENTS C REQUIRED)
find_package (METIS REQUIRED)
find_package(OpenMP)
find_package(OmpSs)

## ASAN is disabled by default enable if needed but
## incompatible with OmpSs (use with care)
set (DISABLE_ASAN TRUE)
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

if (HAVE_FLAG_SANITIZE_ADDRESS AND NOT DISABLE_ASAN)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ASAN_FLAGS} ${OpenMP_C_FLAGS} -O2")
else()
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -O2")
endif()

add_library(backend STATIC
   backend.cpp)

target_include_directories (backend PUBLIC
   ${HDF5_INCLUDE_DIR}
   ${METIS_INCLUDE_DIR})

target_link_libraries (backend PUBLIC
   ${HDF5_LIBRARIES}
   ${METIS_LIBRARY}
   ${CMAKE_DL_LIBS})

add_executable(main_openmp
   main.cpp)

target_link_libraries (main_openmp PUBLIC
   backend)

if(OMPSS_FOUND AND DISABLE_ASAN)

  set (OMPSS_MAIN_SOURCE "${PROJECT_SOURCE_DIR}/ompss_main.cpp")
  set (OMPSS_KERNEL_SOURCE "${PROJECT_SOURCE_DIR}/kernel_ompss.h")

  add_custom_target( ompss ALL
     DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/ompss_main"
  )

  ## OMPSS
  add_custom_command(
     OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/ompss_main"
     COMMAND
       ${OMPSS_CXX} -I${PROJECT_SOURCE_DIR} ${OMPSS_MAIN_SOURCE}
       -o ${CMAKE_CURRENT_BINARY_DIR}/ompss_main
       -L${CMAKE_CURRENT_BINARY_DIR} -lbackend ${METIS_LIBRARY}
       ${HDF5_LIBRARIES} -l${CMAKE_DL_LIBS}
     DEPENDS
       ${OMPSS_MAIN_SOURCE} ${OMPSS_KERNEL_SOURCE} backend
     WORKING_DIRECTORY
       ${PROJECT_SOURCE_DIR}
  )

endif ()

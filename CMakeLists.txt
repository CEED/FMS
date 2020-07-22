# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
# a collaborative effort of two U.S. Department of Energy organizations (Office
# of Science and the National Nuclear Security Administration) responsible for
# the planning and preparation of a capable exascale ecosystem, including
# software, applications, hardware, advanced system engineering and early
# testbed platforms, in support of the nation's exascale computing imperative.

# We use the C_STANDARD target property introduced in CMake 3.1
cmake_minimum_required (VERSION 3.1)

# Project name
project(FMS VERSION 0.1 LANGUAGES C CXX)

# Optional demo - see examples/README.md for more information
option(ENABLE_DEMO "Enable the PUMI + MFEM demo" OFF)
option(FMS_ENABLE_TESTS "Enable FMS testing" OFF)

# Set a default configuration type
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Setup thirdparty dependencies
# Optional conduit implementation for IO functions
if(CONDUIT_DIR)
  find_package(Conduit 0.5.1 REQUIRED
    HINTS ${CONDUIT_DIR}
    NO_DEFAULT_PATH)
  set(FMS_HAVE_CONDUIT TRUE)
  message(STATUS "Using Conduit: ${Conduit_DIR}")
endif()

# For development of the fms-mfem library, will move eventually (still required by demo)
if(MFEM_DIR)
  find_package(mfem REQUIRED
    NAMES MFEM
    HINTS ${MFEM_DIR}
    NO_DEFAULT_PATH)
  set(FMS_HAVE_MFEM TRUE)
  message(STATUS "Using MFEM: ${mfem_DIR}")
endif()

# The following are required for the example
if(ENABLE_DEMO)
  if(CMAKE_VERSION VERSION_LESS 3.10)
    message(FATAL_ERROR "Building the demo requires CMake version 3.10 or later, please rerun cmake with ENABLE_DEMO=OFF")
  endif()

  find_package(MPI QUIET)
  if(MPI_FOUND)
    add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)
    # message(STATUS "Example using MPI_C: ${MPI_C_COMPILER}")
    message(STATUS "Demo using MPI_CXX: ${MPI_CXX_COMPILER}")
    set(FMS_HAVE_MPI TRUE)
  endif()

  if(PUMI_DIR)
    find_package(SCOREC 2.1.0 REQUIRED 
      CONFIG 
      HINTS ${PUMI_DIR} 
      NO_DEFAULT_PATH)
    set(FMS_HAVE_PUMI TRUE)
    message(STATUS "Demo PUMI: ${SCOREC_DIR}")
  endif()
endif()

# Provides cache variables for install dirs
include(GNUInstallDirs)

# Source FMS directory
add_subdirectory(src)

# Examples subdirectory
add_subdirectory(examples)

# Install fmsConfig.cmake
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/fmsConfig.cmake.in
  ${CMAKE_BINARY_DIR}/fmsConfig.cmake
  @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/fmsConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fms)

# Tests/unit tests.
if(FMS_ENABLE_TESTS)
  enable_testing()
  set(INSTALL_GTEST OFF CACHE BOOL "Turn off GTest install")
  add_subdirectory(thirdparty/googletest)
  add_subdirectory(tests)
endif()

# Only enable CPack if we are the top level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(CPACK_PACKAGE_NAME FMS)
  # set(CPACK_PACKAGE_VENDOR MyCompany)
  # set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CPack example project")
  set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
  set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
  set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
  set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
  set(CPACK_VERBATIM_VARIABLES YES)
  set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_LIST_DIR}/README.md)
  set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_LIST_DIR}/LICENSE)
  include(CPack)
endif()

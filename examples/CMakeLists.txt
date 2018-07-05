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

if (ENABLE_DEMO)
  if (NOT PUMI_DIR OR NOT MFEM_DIR)
    message(FATAL_ERROR
            "PUMI_DIR and MFEM_DIR must be specified with ENABLE_DEMO=ON")
  endif()
  # Set MPICXX directly to use an MPI-wrapper, e.g. -DMPICXX=mpicxx
  if (MPICXX)
    set(CMAKE_CXX_COMPILER ${MPICXX})
    # Setting the variables below circumvents autodetection, see FindMPI.cmake.
    set(MPI_CXX_INCLUDE_PATH "")
    set(MPI_CXX_LIBRARIES "")
  endif()
  enable_language(CXX)
  find_package(MPI REQUIRED)
  find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${PUMI_DIR} NO_DEFAULT_PATH)
  find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
               "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)

  add_executable(demo_pumi_mfem demo_pumi_mfem.cpp)
  target_include_directories(demo_pumi_mfem PUBLIC ${MPI_CXX_INCLUDE_PATH}
                             ${MFEM_INCLUDE_DIRS})
  target_link_libraries(demo_pumi_mfem fms SCOREC::core mfem
                        ${MPI_CXX_LIBRARIES})

  # Fetch and/or copy the data directory, if necessary
  add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/data
    COMMAND git clone https://github.com/mfem/data.git data
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    COMMENT "Fetching MFEM-data from GitHub to the source tree ...")
  add_custom_target(fetch_data DEPENDS ${PROJECT_SOURCE_DIR}/data)
  add_dependencies(demo_pumi_mfem fetch_data)
  if (NOT ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}"))
    add_custom_command(OUTPUT data_is_copied
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              ${PROJECT_SOURCE_DIR}/data ${PROJECT_BINARY_DIR}/data
      COMMAND ${CMAKE_COMMAND} -E touch data_is_copied
      COMMENT "Copying the data directory to the build tree ...")
    add_custom_target(copy_data DEPENDS fetch_data data_is_copied)
    add_dependencies(demo_pumi_mfem copy_data)
  endif()
endif()

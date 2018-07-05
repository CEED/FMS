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
project(FMS VERSION 0.1 LANGUAGES C)

# Optional demo - see examples/README.md for more information
option(ENABLE_DEMO "Enable the PUMI + MFEM demo" OFF)

# Set a default configuration type
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Source FMS directory
add_subdirectory(src)

# PUMI_DIR and MFEM_DIR are needed when ENABLE_DEMO=ON
set(PUMI_DIR "" CACHE PATH "Path to the PUMI library.")
set(MFEM_DIR "" CACHE PATH "Path to the MFEM library.")

# Examples subdirectory
add_subdirectory(examples)

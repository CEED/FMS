# Copyright (c) 2020, Lawrence Livermore National Security, LLC. Produced at
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

set(FMS_SOURCES fms.c fmsio.c)
set(FMS_HEADERS fms.h fmsio.h)

add_library(fms ${FMS_SOURCES} ${FMS_HEADERS})
set_target_properties(fms
  PROPERTIES
    C_STANDARD 99
    C_EXTENSIONS OFF)

target_include_directories(fms PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)

# The IO functions will be optionally implemented with Conduit
if(FMS_HAVE_CONDUIT)
  target_link_libraries(fms PRIVATE conduit::conduit)
  target_compile_definitions(fms PRIVATE FMS_HAVE_CONDUIT=1)
  set_target_properties(fms PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# Install main fms library
install(TARGETS fms DESTINATION ${CMAKE_INSTALL_LIBDIR} EXPORT fms)
install(FILES ${FMS_HEADERS} DESTINATION include)
install(EXPORT fms 
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fms
  NAMESPACE fms::
  FILE fms.cmake)

# Setup fms-mfem library which ultimately may not live in this repo
# if(FMS_HAVE_MFEM)
#   add_library(fms-mfem fms-mfem.cpp fms-mfem.h)
#   set_property(TARGET fms-mfem PROPERTY CXX_STANDARD 11)
#   target_include_directories(fms-mfem 
#     PUBLIC
#       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#       $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include)
#   target_link_libraries(fms-mfem PUBLIC mfem)

#   install(TARGETS fms-mfem DESTINATION ${CMAKE_INSTALL_LIBDIR} EXPORT fms-mfem)
#   install(FILES fms-mfem.h DESTINATION include)
#   install(EXPORT fms-fmem 
#     DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fms
#     NAMESPACE fms::
#     FILE fms-mfem.cmake)
# endif()
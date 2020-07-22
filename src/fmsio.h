/*
 Copyright (c) 2020, Lawrence Livermore National Security, LLC. Produced at
 the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
 reserved. See files LICENSE and NOTICE for details.

 This file is part of CEED, a collection of benchmarks, miniapps, software
 libraries and APIs for efficient high-order finite element and spectral
 element discretizations for exascale applications. For more information and
 source code availability see http://github.com/ceed.

 The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
 a collaborative effort of two U.S. Department of Energy organizations (Office
 of Science and the National Nuclear Security Administration) responsible for
 the planning and preparation of a capable exascale ecosystem, including
 software, applications, hardware, advanced system engineering and early
 testbed platforms, in support of the nation's exascale computing imperative.
*/

#ifndef FMS_IO_HEADER
#define FMS_IO_HEADER
#include <fms.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
@brief Writes the provided FmsDataCollection to a file.
@param filename The name of the file to save. If the filename does not include 
                an appropriate file extension, one may be added.
@param protocol The type of file to be saved. By default, the protocol
                should be "ascii". If FMS is compiled with Conduit support then
                the protocol can match the supported Conduit protocols, which
                can include: "json", "yaml", "silo", "hdf5".
@param dc       The FMS object that will be written to the file.
@return The function returns 0 on success and non-zero for failure.
*/
int FmsIOWrite(const char *filename, const char *protocol, FmsDataCollection dc);

/**
@brief Reads an FmsDataCollection from a file.
@param filename The name of the file to read.
@param protocol The type of file to be read. By default, the protocol
                should be "ascii". If FMS is compiled with Conduit support then
                the protocol can match the supported Conduit protocols, which
                can include: "json", "yaml", "silo", "hdf5". Adding a protocol
                will help FMS and Conduit decide which method is appropriate for
                reading the file.
@param[out] dc  The FMS object that was read from the file.
@return The function returns 0 on success and non-zero for failure.
*/
int FmsIORead(const char *filename, const char *protocol, FmsDataCollection *dc);


#ifdef __cplusplus
} // extern "C"
#endif

#endif // FMS_IO_HEADER
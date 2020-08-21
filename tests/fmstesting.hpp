#ifndef FMS_TESTING_HPP
#define FMS_TESTING_HPP

#include <fms.h>

// 1 element quad mesh
int Construct2DData0(FmsDataCollection *out_dc);

// 1 element triangle mesh
int Construct2DData1(FmsDataCollection *out_dc);

// 1 element hex mesh
int Construct3DData0(FmsDataCollection *out_dc);

// 1 element tet mesh
int Construct3DData1(FmsDataCollection *out_dc);

#endif

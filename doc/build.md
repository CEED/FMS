# Building FMS

FMS requires CMake (atleast version 3.1) and a C99 compiler to build.

Optionally Conduit can be linked in to use conduit_relay for FmsIO operations. 
This enables FmsIO to write conduit binary or other conduit_relay supported protocols.

## Default configuration

FMS will build optimized for release and as a static library if you use the following command line

```sh
mkdir build
cd build
cmake ..
make -j 2
```

## Example - building in debug mode as a shared library with testing enabled

```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_INSTALL_PREFIX=../install -DFMS_ENABLE_TESTS=ON ..
make -j 4
ctest -V # or ctest --output-on-failure
make install
```

## CMake options

|Option|Description|Default Value|
|------|-----------|-------------|
|CMAKE_BUILD_TYPE|Can be set to Debug, Release, RelWithDebInfo, or MinSizeRel.<br>Example: *-DCMAKE_BUILD_TYPE:STRING=Release*|Release|
|CMAKE_INSTALL_PREFIX|The desired install location to be used by "make install".<br>Example: *-DCMAKE_INSTALL_PREFIX:PATH=/usr/local/fms*|Not set|
|BUILD_SHARED_LIBS|Whether or not to build Fms as a shared library.<br>Example: *-DBUILD_SHARED_LIBS=ON*|Not set|
|CONDUIT_DIR|Prefix of your Conduit installation.<br>Example: *-DCONDUIT_DIR=/usr/local/conduit*|Not set|
|FMS_ENABLE_TESTS|Enables unit tests. This requires a C++11 compiler and will add googletest to the build.|OFF|
|FMS_ENABLE_DEMO|Enables the pumi/mfem demo located in examples/demo_pumi_mfem.<br>Enabling this demo **requires** that you specify PUMI_DIR and MFEM_DIR|OFF|
|PUMI_DIR|Prefix of your pumi/scorec-core installation.<br>Example: *-DPUMI_DIR=/usr/local/SCOREC*|Not set|
|MFEM_DIR|Prefix of your mfem installation.<br>Example: *-DMFEM_DIR=/usr/loca/mfem*|Not set|

*Note*: Checkout the demo_pumi_mfem [README](../examples/demo_pumi_mfem/README.md) for more info on building/running the demo.

## Linking to Fms

When installed Fms will create all of the files needed to link to it with CMake.
There is an example of how todo this in examples/include-fms.
After installing Fms you can run the following command to test this:

```sh
cd examples/include-fms
mkdir build
cd build
cmake -DFMS_DIR=*fms/install/prefix* ..
make
./main
```
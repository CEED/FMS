## FMS Examples

### PUMI + MFEM demo

The demo performs the following steps:
* loads a PUMI mesh
* converts the PUMI mesh object to an FMS mesh object
* converts the FMS mesh object to an MFEM mesh object
* sends the MFEM mesh object to GLVis for visualization.

To build the demo, first clone and build the PUMI and MFEM repositories, e.g. next to the FMS directory:
```console
git clone https://github.com/SCOREC/core.git pumi
git clone https://github.com/mfem/mfem.git
cd pumi
mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
  -DSCOREC_CXX_WARNINGS=OFF -DCMAKE_INSTALL_PREFIX=../../pumi-install
make -j 4 install
cd ../..
cd mfem
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../mfem-install
make -j 4 install
cd ../..
```
Next, fetch the MFEM data repository used by the demo:
```console
cd fms
git clone https://github.com/mfem/data.git
```
And finally, build FMS and the demo:
```console
mkdir build
cd build
cmake .. -DENABLE_DEMO=ON -DMPICXX=mpicxx \
  -DPUMI_DIR=../../pumi-install -DMFEM_DIR=../../mfem-install
make
cd examples
./demo_pumi_mfem
```
The last step in the demo will send the mesh to [GLVis](http://glvis.org/) for visualization, if a GLVis server is running. GLVis can be cloned and built using the following commands, e.g. next to the FMS, PUMI, and MFEM directories (see the [glvis/INSTALL](https://github.com/GLVis/glvis/blob/master/INSTALL) file for additional required libraries):
```console
git clone https://github.com/GLVis/glvis.git
cd glvis
mkdir build
cd build
cmake .. -DMFEM_DIR=../../mfem-install
make -j 4
```
To start GLVis in server mode use:
```console
./glvis
```

ARC RICH detector V0
=======================

Simple cell description of the ARC RICH detector.

# Directory structure

README.md: this file
CMakeLists.txt: minimal cmake macro to configure the detector construction
compact: has material and elements with their description (including optical properties). The compact file of the detector is also here.
src: detector consctructor, writen in c++
myscripts: useful bash scripts and GEANT4 macro files


# Basic commands

Build with `cmake`, 
```bash
source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
cmake -B build -S . -D CMAKE_INSTALL_PREFIX=install
cmake --build build -- install
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
```

To display the geometry

```bash
geoDisplay compact/arc_v0.xml 
```

To show materials in along a given line (in this case from origin to (0,0,-100)cm)

```bash
materialScan compact/arc_v0.xml 0 0 0 0 0 -100
```

Compile and install after every modification of the c++ detector constructor code.

```bash
cmake --build build -- install
```

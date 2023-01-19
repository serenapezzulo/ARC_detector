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

To check overlaps, we can use Geant4 check
```shell
ddsim --compactFile ./compact/arc_v0.xml --runType run --part.userParticleHandler='' --macroFile myscripts/overlap.mac
```

# Run simple simulation

Use this command to open the Geant4 Qt application,

```shell
ddsim --compactFile compact/arc_v0.xml \ 
      --enableGun \ 
      --gun.distribution uniform \ 
      --gun.energy "10*GeV" \ 
      --gun.particle mu- \ 
      --numberOfEvents 1 \ 
      --runType qt \ 
      --part.userParticleHandler='' \ 
      --macroFile myscripts/vis.mac 
```

then press play button on top to launch one event.

# Useful links

Documentation of DD4hep,
* PDF: https://dd4hep.web.cern.ch/dd4hep/usermanuals/DD4hepManual/DD4hepManual.pdf
* Code: https://github.com/AIDASoft/DD4hep
* Class reference: https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html
* Beginners guide: https://dd4hep.web.cern.ch/dd4hep/page/beginners-guide/
* FCC tutorial: https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Geometry/Geometry.html#geometry-driver-modifications
* NPdet: https://eic.phy.anl.gov/npdet/index.html

Example of RICH detector implemented in DD4hep
* Proximity Focusing RICH: https://github.com/AIDASoft/DD4hep/tree/master/examples/OpticalTracker

Documentation of ARC,
* Dec-22 ARC meeting indico: https://indico.cern.ch/event/1231098/
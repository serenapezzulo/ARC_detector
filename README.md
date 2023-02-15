ARC RICH detector V0
=======================

Simple cell description of the ARC RICH detector.

The directory structure is the following

* README.md: this file
* CMakeLists.txt: minimal cmake macro to configure the detector construction
* compact: has material and elements with their description (including optical properties). The compact file of the detector is also here.
* src: detector consctructor, writen in c++
* myscripts: useful bash scripts and GEANT4 macro files


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

To convert the geometry into ROOT

```bash
./myscripts/dd4hep2root -c  compact/arc_v0.xml -o arc_v0.root
```

# Things to do after every change

Compile and install after every modification of the c++ detector constructor code.

```bash
cmake --build build -- install
```

To show materials in along a given line (in this case from origin to (0,0,-100)cm), check if the modified volume has the proper material

```bash
materialScan compact/arc_v0.xml 0 0 0 0 0 -100
```

Check overlaps, using Geant4 check
```shell
ddsim --compactFile ./compact/arc_v0.xml --runType run --part.userParticleHandler='' --macroFile myscripts/overlap.mac >> overlapDump.txt
```


# Run simulation

The script `arcsim.py` is a ddsim steering file, which enable Qt visualization, setup the cerenkov physic lists and the default input and output files. The configuration inside the steering file can be overriden if options are added after the name of the steering file (see -N 1 example below). 

Use this command to open the Geant4 Qt application,

```shell
python arcsim.py -N 1
```

then press play button on top to launch one event. The geometry and tracks of one charged pion (pink) crossing the detector, and the cerenkov photons (grey) are displayed below.


![Visualization of geometry and tracks of one charged pion (pink) crossing the detector, and the cerenkov photons (grey).](https://mattermost.web.cern.ch/files/4eyjw6q8b38tfkmh5nrf7qytsr/public?h=rLXJPn3Zsd6-0g3Q9ZaqHI9pAFglRqD_kix70QP6nXs)



More events can be simulated, but in this case only 1 will be saved into the output file `arcsim.root`. To visualize the hit pattern in the detector, we can open the root file and write the following in the command line

```cpp
EVENT->Draw("ARC_HITS.position.Y():ARC_HITS.position.X()");
```

A new canvas will open showing something similar to the following graph:

![Hit pattern of photons (and pion) in the detector](https://mattermost.web.cern.ch/files/11f17b5nctfkjqjsw1cmoh1rko/public?h=OQCOs1RkwC560pOU0reOnG9pJabN4rDqTu2wgqHeHNg)

Hit pattern of photons (and pion) in the detector


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
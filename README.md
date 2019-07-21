# [TURTLE][TURTLE]-perfs [![Build Status](https://travis-ci.com/niess/turtle-perfs.svg?branch=master)](https://travis-ci.com/niess/turtle-perfs)
_Performance tests of the [TURTLE][TURTLE] library._

This repository contains a collection of utilities and scripts for testing the
performances of the [TURTLE][TURTLE] library when navigating through a
topography represented by a Digital Elevation Model (DEM). The detailed results
of these tests can be found in [arXiv:1904.03435][ARXIV]. In addition to the
_optimistic_ algorithm included in [TURTLE][TURTLE] various alternative ray
tracing algorithms are used for purpose of comparison. Based on these
ray tracers four geometry wrappers of topography data have been implemented.
The corresponding source files are:

1. [src/bvh.cpp](src/bvh.cpp): the default BVH accelerated geometry using
   [CGAL][CGAL].

2. [src/embree.c](src/embree.c): an alternative BVH accelerated geometry
   using [Embree][EMBREE].
   > This is 10 times faster than the [CGAL][CGAL] version but much less accurate
   > as well. It uses two times less memory per topography node.

3. [src/mesh.cpp](src/mesh.cpp): a polyhedral mesh using [CGAL](CGAL) for
   geometry operations.

4. [src/opti.c](src/opti.c): the _optimistic_ algorithm using a
   [`turtle_stepper`][TURTLE_STEPPER].

The generic ray tracing algorithm is implemented in
[src/tracer.c](src/tracer.c). The [PUMAS][PUMAS] based simulation test is
located in [src/simulator.c](src/simulator.c). Note that [PUMAS][PUMAS] is
[patched](src/patch/pumas.patch) in order to force all particles to follow
a straight trajectory, even though detailed scattering is simulated.


## Installation

The test suite was designed for Linux. There is no binary distribution. It must
be built and installed from the source using [CMake][CMAKE] 3 or later. Most
external dependencies will be automatically downloaded and installed as well.
Yet the following executables and libraries are assumed to be already installed
on your system:

- [`libboost-system-dev`][BOOST] and [`libboost-thread-dev`][BOOST] for
  [CGAL][CGAL].
- [`libpng-dev`][PNG] for topography maps with [TURTLE][TURTLE].
- [wget][WGET] for downloading the Topography data. You can also do this
  manually from the [releases][RELEASES] page.

Installation can be done as following:
```bash
git clone https://github.com/niess/turtle-perfs.git && cd turtle-perfs
mkdir build && cd build
cmake .. && make install
```
By default the installation is done to the source directory with a `Release`
build.  Building the [Geant4][GEANT4] test suite is optional with the
`USE_GEANT` [CMake][CMAKE] flag. It requires an existing installation of
[Geant4][GEANT4].  This can be done as:
```
source /path/to/geant4/bin/geant4.sh
cmake .. -DUSE_GEANT=ON && make install
```


## Usage

The test environment can be initialised by sourcing the provided
[`setup.sh`](setup.sh) files, as:
```bash
source setup.sh
```
The tests scripts are located under the [`scripts`](scripts) folder. Below is
a listing of the available scripts as well as examples of usage illustrating
their Command Line Interface (CLI). Extra informations can be obtained with the
`-h` option of the CLI.

---
- [trace.lua](scripts/trace.lua): do a ray tracing from a given point of view
  and along a given line of sight. The three algorithms are tested (BVH,
  polyhedral mesh & _optimistic_). For example:
  ```bash
  trace.lua CDC 26 10 --period=10 --embree
  ```
  does a ray tracing test starting from the CDC view and going along an azimuth
  of 26 deg and an elevation of 10 deg, with a map down-sampled by a factor of
  10 (in _x_ and _y_ axis) and with Embree's BVH algorithm (default is CGAL's
  BVH).

---
- [scan-trace.lua](scripts/scan-trace.lua): do a rock depth scan from the
  given point of view and with the given algorithm. The down-sampling period
  of the DEM must be provided as well. For example:
  ```bash
  scan-trace.lua ULS opti 1
  ```
  does a rock depth scan from the Ulastai view point using the _optimistic_
  algorithm and the full DEM.

---
- [simulate.lua](scripts/simulate.lua): run a muography simulation from a given
  point of view and along a given line of sight. The three ray tracing
  algorithms are tested (BVH, polyhedral mesh & _optimistic_). For example:
  ```bash
  simulate.lua ULS 90 5 --period=100 --precompute
  ```
  does the simulation test with a map down-sampled by a factor of 100 along _x_
  and _y_ axes and with a pre-computed geometry.

---
- [scan-simulate.lua](scripts/scan-simulate.lua): run a muography scan from the
  given point of view and with the given ray tracing algorithm. The
  down-sampling period of the DEM must be provided as well. For example:
  ```bash
  scan-simulate.lua CDC BVH 30 --embree
  ```
  does the simulation from the Ulastai view point using Embree's rau tracing
  and a DEM down-sampled by 30 along _x_ and _y_ axes.

---
- [geant4.lua](scripts/geant4.lua): run a [Geant4][GEANT4] simulation from a
  given point of view and along a given line of sight. The default simulation
  generates a 10 TeV `Geantino` primary and uses a [G4Turtle][G4TURTLE]
  geometry. This can be modified by command line options. For example:
  ```bash
  geant4.lua CDC 26 10 --particle=MuonMinus --tessellate
  ```
  generates a 10 TeV primary &mu; and uses a G4TessellatedSolid geometry. 

---
- [scan-geant4.lua](scripts/scan-geant4.lua): run a [Geant4][GEANT4] scan from
  the given point of view and with the given down-sampling. For example:
  ```bash
  scan-geant4.lua CDC 1
  ```
  runs a `Geantino` scan with a G4Turtle geometry. As previously, the primary
  and the geometry can be overridden by command line options.


## License

The tests procedures are [unlicensed](LICENSE), i.e. freely available to the
public domain. They can be copied, modified, etc. without any notice. The
external dependencies have their own specific licences.

[ARXIV]: https://arxiv.org/abs/1904.03435
[BOOST]: https://www.boost.org/
[CGAL]: https://www.cgal.org/
[CMAKE]: https://cmake.org/
[EMBREE]: https://www.embree.org/
[GEANT4]: https://geant4.web.cern.ch/
[G4TURTLE]: https://github.com/niess/turtle-geant4
[GEANT4]: https://geant4.web.cern.ch/
[PNG]: http://www.libpng.org/pub/png/libpng.html
[PUMAS]: https://niess.github.io/pumas-pages/
[RELEASES]: https://github.com/niess/turtle-perfs/releases/tag/topography
[TURTLE]: https://niess.github.io/turtle-pages/
[TURTLE_STEPPER]: https://niess.github.io/turtle-docs/#HEAD/group/stepper
[WGET]: https://www.gnu.org/software/wget/

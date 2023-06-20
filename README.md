# NIM-CONV22
The solution of the time-dependent, incompressible, Navier-Stokes equations, coupled with the temperature equation in curvilinear coordinates using nodal integral method. <br />
<br />
The code is written by Ibrahim Jarrah. If you have any questions, feel free to reach out at [ijarrah@anl.gov](mailto:ijarrah@anl.gov).
# Test case: Natural Convection in a Cylindrical Enclosure with an Internal Hexagonal Object
The details of the case can be found in this  [paper](https://www.dl.begellhouse.com/references/1bb331655c289a0a,715716a453fb6732,6df126bb2596ded2.html).<br />
# Silo Library
[Silo](https://github.com/LLNL/Silo) library is required to run the code. You can istall it and provide the path to silo in the Makefile, or use the provided script [install.sh](3rd_part/install.sh).
# Installation
```
git clone https://github.com/ijarrah2/NIM-CONV22.git
cd NIM-CONV22/3rd_party
./install.sh
cd ..
```
# Running the case
```
./run.sh
```
# Running the case with different parameters
To change the test case setup, modify the [src/parameters.f90](src/parameters.f90) file.
# Visualization
Geometry and results can be visualized using [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit). The code create binary silo outputs with a master file (plots.visit) that can be opened using VisIt.

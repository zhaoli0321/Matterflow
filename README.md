# Matterflow: README
## Introduction

The Lagrangian hydrodynamics are commonly used to simulate problems such as detonation, inertial confinement fusion, hypervelocity collision. A classical Lagrangian staggered-grid hydrodynamic (SGH) method easily yields the cell-to-cell oscillations of physical quantities on triangular grids  (denoted as ``checkerboard oscillations") . A novel flux method, called the matter flow method, is proposed to alleviate the checkerboard oscillations for the SGH method. The method is integrated in the *Matterflow* software package.

*Matterflow* is based on the flux method described in the following article:

> Li Zhao, Bo Xiao, Ganghua Wang, Haibo Zhao, Jinsong Bai, Chunsheng Feng, Shi Shu,
>
> "A Matter Flow Method for Staggered Lagrangian Hydrodynamics on Triangular Grids".

## Directory Structure

- data : This folder contains data (meshes) files;
- include : This folder contains header files;
- lib : This folder contains library file;
- main : This folder contains main function of the Matterflow method;
- src  :  This folder contains sources code;
- util :  This folder contains tools - automatically generate header files.

## Build

Build *Matterflow*:

```makefile
~> tar -zxvf Matterflow.tar.gz
~> cd Matterflow
~> make 
```

## Running

Taylor-Green vortex grid files (Unstructured grids: h=0.05, h=0.025, h=0.0125) are located in the data directory.

```sh
./test.ex -f ./data/0.05/Compute_Model_For_TriAngels.txt
./test.ex -k 2 -cs 0.5 -f ./data/0.025/Compute_Model_For_TriAngels.txt
./test.ex -help # Enter -help/--help/-h to output optionals.
```

## Results

The terminal outputs:

```sh
Loading compute model information from file "./data/0.05/Compute_Model_For_TriAngels.txt" 

Compute model information loaded. 

HAVE_MF = 1, TotBitmapOutputs = 20 

Times of viscosity coefficient 2.0000, C_safe = 0.5000

Domain Omega = [0, 1.00] x [0, 1.00]

MinMeshScale = 0.034082, MaxMeshScale = 0.063082, AvgMeshScale = 0.048582 

Evolving: 
  0: Iter =      0, t = 0.000000
  1: Iter =     36, t = 2.537650e-02, deltT = 7.025646e-04, TimeUsed=0.09s
  2: Iter =     72, t = 5.054822e-02, deltT = 6.957280e-04, TimeUsed=0.16s
  3: Iter =    112, t = 7.547126e-02, deltT = 5.587029e-04, TimeUsed=0.22s
  4: Iter =    162, t = 1.002846e-01, deltT = 4.447163e-04, TimeUsed=0.30s
  5: Iter =    225, t = 1.251652e-01, deltT = 3.535268e-04, TimeUsed=0.40s
  6: Iter =    305, t = 1.502629e-01, deltT = 2.804681e-04, TimeUsed=0.52s
  7: Iter =    405, t = 1.751864e-01, deltT = 2.229943e-04, TimeUsed=0.67s
  8: Iter =    531, t = 2.001671e-01, deltT = 1.773990e-04, TimeUsed=0.87s
  9: Iter =    689, t = 2.251146e-01, deltT = 1.413826e-04, TimeUsed=1.10s
 10: Iter =    887, t = 2.500580e-01, deltT = 1.128899e-04, TimeUsed=1.39s
 11: Iter =   1135, t = 2.750306e-01, deltT = 9.030808e-05, TimeUsed=1.77s
 12: Iter =   1445, t = 3.000319e-01, deltT = 7.240004e-05, TimeUsed=2.25s
 13: Iter =   1831, t = 3.250265e-01, deltT = 5.820429e-05, TimeUsed=2.82s
 14: Iter =   2310, t = 3.500034e-01, deltT = 4.694058e-05, TimeUsed=3.52s
 15: Iter =   2904, t = 3.750220e-01, deltT = 3.797020e-05, TimeUsed=4.38s
 16: Iter =   3636, t = 4.000124e-01, deltT = 3.083575e-05, TimeUsed=5.44s
 17: Iter =   4536, t = 4.250166e-01, deltT = 2.514227e-05, TimeUsed=6.76s
 18: Iter =   5637, t = 4.500162e-01, deltT = 2.059465e-05, TimeUsed=8.36s
 19: Iter =   6977, t = 4.750050e-01, deltT = 1.695677e-05, TimeUsed=10.31s
 20: Iter =   8601, t = 5.000103e-01, deltT = 1.403848e-05, TimeUsed=12.66s


                            Numerical Simulation Information                                

--------------------------------------------------------------------------------------------
The number of triangles: 1044, the number of vertices: 563, Iter = 8601
Simulation Time:                 12.6600
MatterFlow Time:                 8.8000 (69.51%)
--------------------------------------------------------------------------------------------
```


The test results are saved in the "results" directory:

- TecPlot: The TecPlot files are used to display contour plots;
- Physical_Quantity_Statistics.txt: Physical quantity statistics;
- etc.

## License

This software is free software distributed under the Lesser General Public License or LGPL, version 3.0 or any later versions. This software distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program. If not, see http://www.gnu.org/licenses/.


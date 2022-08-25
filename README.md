# Matterflow: README

## Introduction

The Lagrangian hydrodynamics are commonly used to simulate problems such as detonation, inertial confinement fusion, hypervelocity collision. A classical Lagrangian staggered-grid hydrodynamic (SGH) method easily yields the cell-to-cell oscillations on triangular grids. A novel flux method, so-called the matter-flow method, is proposed to eliminate such oscillations for the SGH method. The proposed method is implemented in the *Matterflow* software package.

*Matterflow* is based on the flux method described in the following article:

> Li Zhao, Bo Xiao, Ganghua Wang, Haibo Zhao, Jinsong Bai, Chunsheng Feng, Shi Shu,
>
> "Study on a Matter Flux Method for Staggered Essentially Lagrangian Hydrodynamics on Triangular Grids".

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

The grid files are located in the data directory.

```sh
for example:
./test.ex -tg 1 -k 0.5 -mf 1 -f ./data/Taylor-Green/0.05/Compute_Model_For_TriAngels.txt
./test.ex -help # Enter -help or --help or -h to output optionals.
```

## Results

The terminal outputs:

```sh
Loading compute model information from file "./data/Taylor-Green/0.05/Compute_Model_For_TriAngels.txt" 

Compute model information loaded. 

HAVE_MF = 1, TotBitmapOutputs = 20 

Times of viscosity coefficient 0.5000

Times of heat diffusion coefficient 0.0000

Time-step factor C_safe = 0.1000

Domain Omega = [0, 1.00] x [0, 1.00]

MinMeshScale = 0.034082, MaxMeshScale = 0.063082, AvgMeshScale = 0.048582 

Evolving: 
  0: Iter =      0, t = 0.000000
  1: Iter =     43, t = 2.558046e-02, deltT = 5.935042e-04, TimeUsed=0.36s
  2: Iter =     85, t = 5.043139e-02, deltT = 5.899141e-04, TimeUsed=0.43s
  3: Iter =    130, t = 7.518717e-02, deltT = 5.116053e-04, TimeUsed=0.49s
  4: Iter =    183, t = 1.002545e-01, deltT = 4.393897e-04, TimeUsed=0.56s
  5: Iter =    244, t = 1.250702e-01, deltT = 3.783827e-04, TimeUsed=0.64s
  6: Iter =    316, t = 1.502880e-01, deltT = 3.255884e-04, TimeUsed=0.73s
  7: Iter =    398, t = 1.750626e-01, deltT = 2.814595e-04, TimeUsed=0.84s
  8: Iter =    494, t = 2.001410e-01, deltT = 2.430743e-04, TimeUsed=0.96s
  9: Iter =    605, t = 2.251272e-01, deltT = 2.091574e-04, TimeUsed=1.10s
 10: Iter =    734, t = 2.501246e-01, deltT = 1.800978e-04, TimeUsed=1.26s
 11: Iter =    883, t = 2.750074e-01, deltT = 1.553118e-04, TimeUsed=1.46s
 12: Iter =   1057, t = 3.000663e-01, deltT = 1.339257e-04, TimeUsed=1.67s
 13: Iter =   1258, t = 3.250481e-01, deltT = 1.156535e-04, TimeUsed=1.92s
 14: Iter =   1491, t = 3.500709e-01, deltT = 9.998566e-05, TimeUsed=2.22s
 15: Iter =   1760, t = 3.750694e-01, deltT = 8.659641e-05, TimeUsed=2.56s
 16: Iter =   2070, t = 4.000496e-01, deltT = 7.517478e-05, TimeUsed=2.95s
 17: Iter =   2427, t = 4.250554e-01, deltT = 6.542578e-05, TimeUsed=3.39s
 18: Iter =   2836, t = 4.500301e-01, deltT = 5.712965e-05, TimeUsed=3.92s
 19: Iter =   3304, t = 4.750286e-01, deltT = 5.006250e-05, TimeUsed=4.49s
 20: Iter =   3837, t = 5.000258e-01, deltT = 4.403802e-05, TimeUsed=5.15s


                            Numerical Simulation Information                              
--------------------------------------------------------------------------------------------
The number of triangles: 1044, the number of vertices: 563, Iter = 3837
Simulation Time:                 5.1600
MatterFlow Time:                 3.2000 (62.02%)
--------------------------------------------------------------------------------------------

Print TaylorGreenL2Norm
Density  L2 = 1.143512e-03
Pressure L2 = 1.694569e-02
Velocity L2 = 9.235076e-03
```

The test results are saved in the "results" directory:

- TecPlot: The TecPlot files are used to display contour plots;
- Physical_Quantity_Statistics.txt: Physical quantity statistics;
- etc.

## License

This software is free software distributed under the Lesser General Public License or LGPL, version 3.0 or any later versions. This software distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program. If not, see http://www.gnu.org/licenses/.

# Matterflow: README

## Introduction

The Lagrangian hydrodynamics are commonly used to simulate problems such as detonation, inertial confinement fusion, hypervelocity collision. A classical Lagrangian staggered-grid hydrodynamic (SGH) method easily yields the cell-to-cell oscillations on triangular grids. A novel flux method, so-called the matter-flow method, is proposed to eliminate such oscillations for the SGH method. The proposed method is implemented in the *Matterflow* software package.

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

The grid files are located in the data directory.

```sh
for example:
./test.ex -tg 1 -k 1 -cs 0.5 -mf 1 -bn 20 -f ./data/Taylor-Green/0.05/Compute_Model_For_TriAngels.txt
./test.ex -help # Enter -help or --help or -h to output optionals.
```

## Results

The terminal outputs:

```sh
Loading compute model information from file "./data/Taylor-Green/0.05/Compute_Model_For_TriAngels.txt" 

Compute model information loaded. 

HAVE_MF = 1, TotBitmapOutputs = 20 

Times of viscosity coefficient 1.0000

Times of heat diffusion coefficient 0.0000

Time-step factor C_safe = 0.5000

Domain Omega = [0, 1.00] x [0, 1.00]

MinMeshScale = 0.034082, MaxMeshScale = 0.063082, AvgMeshScale = 0.048582 

Evolving: 
  0: Iter =      0, t = 0.000000
  1: Iter =      9, t = 2.677161e-02, deltT = 2.966986e-03, TimeUsed=0.29s
  2: Iter =     17, t = 5.044303e-02, deltT = 2.952731e-03, TimeUsed=0.31s
  3: Iter =     26, t = 7.503122e-02, deltT = 2.556436e-03, TimeUsed=0.33s
  4: Iter =     37, t = 1.007926e-01, deltT = 2.181343e-03, TimeUsed=0.35s
  5: Iter =     50, t = 1.268001e-01, deltT = 1.861921e-03, TimeUsed=0.37s
  6: Iter =     64, t = 1.508989e-01, deltT = 1.611252e-03, TimeUsed=0.40s
  7: Iter =     81, t = 1.761829e-01, deltT = 1.388204e-03, TimeUsed=0.43s
  8: Iter =    100, t = 2.006417e-01, deltT = 1.205171e-03, TimeUsed=0.46s
  9: Iter =    122, t = 2.252301e-01, deltT = 1.043844e-03, TimeUsed=0.50s
 10: Iter =    148, t = 2.503404e-01, deltT = 9.004991e-04, TimeUsed=0.54s
 11: Iter =    178, t = 2.753701e-01, deltT = 7.783454e-04, TimeUsed=0.59s
 12: Iter =    213, t = 3.006225e-01, deltT = 6.729622e-04, TimeUsed=0.64s
 13: Iter =    252, t = 3.250405e-01, deltT = 5.857172e-04, TimeUsed=0.70s
 14: Iter =    298, t = 3.500962e-01, deltT = 5.091388e-04, TimeUsed=0.77s
 15: Iter =    351, t = 3.752272e-01, deltT = 4.436587e-04, TimeUsed=0.84s
 16: Iter =    411, t = 4.000810e-01, deltT = 3.881851e-04, TimeUsed=0.93s
 17: Iter =    481, t = 4.253269e-01, deltT = 3.366387e-04, TimeUsed=1.03s
 18: Iter =    560, t = 4.501458e-01, deltT = 2.944729e-04, TimeUsed=1.15s
 19: Iter =    651, t = 4.752393e-01, deltT = 2.593356e-04, TimeUsed=1.29s
 20: Iter =    753, t = 5.001532e-01, deltT = 2.310233e-04, TimeUsed=1.43s


                            Numerical Simulation Information                              
--------------------------------------------------------------------------------------------
The number of triangles: 1044, the number of vertices: 563, Iter = 753
Simulation Time:                 1.4400
MatterFlow Time:                 0.6800 (47.22%)
--------------------------------------------------------------------------------------------

Print TaylorGreenL2Norm
Density  L2 = 1.415778e-03
Pressure L2 = 1.428664e-02
Velocity L2 = 1.159701e-02
```

The test results are saved in the "results" directory:

- TecPlot: The TecPlot files are used to display contour plots;
- Physical_Quantity_Statistics.txt: Physical quantity statistics;
- etc.

## License

This software is free software distributed under the Lesser General Public License or LGPL, version 3.0 or any later versions. This software distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program. If not, see http://www.gnu.org/licenses/.

/*! \file  GlobalParameters.c
 *
 *  \brief Define some global parameters
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"


/// Model file name
char *modelInfoFileName = "./data/0.05/Compute_Model_For_TriAngels.txt";

/// Whether to use the matter flow method, 1: ON, 0: OFF
int HAVE_MF     = 1;     

/// Whether it is a Taylor-Green vortex Problem, the energy equation needs to add a source term
int TaylorGreen = 1;     

/// Times of viscosity coefficient (adjustable)
double ViscCoeffTimes = 2; 

/// Time-step factor
double timeStepSafeFactor = 0.5; 

/// Gas Universal Constant
double ROfGas = 8.3149 * 1.0e-5; // 8.3149*J.mol^-1.k^-1

/// Records the time of system evolution for continued calculation.
double TimeEvolvedRecord = 0.0;

/// The number of times the system has evolved. (it is equal to time divided by some interval)
int TimeNumbers = 0;

/// Number of pictures (Tecplot files) (it is equal to time divided by picture output interval)
int TimeCountForBitmap = 0;

/// Time interval used to control data output
double TimeInterval = 1.0;

/// Used to control the end time of the program running
double TimeEnd = 1.0;

/// The number of pictures (Tecplot files) to output
int TotBitmapOutputs = 20;

/// Initialize the number of iterations
int Iter = 0;

/// The file name of the statistical physical quantites
FILE * PhysQuantStatisticsFile;

/// Total simulation time (initialization)
double TotTime = 0; 

/// Matter flow time (initialization)
double MatterFlowTime = 0;

/// Auxiliary variable (initialization)
clock_t Start_time = 0;

/// Auxiliary variable (initialization)
clock_t End_time = 0;

/// Auxiliary variable (initialization)
bool CellDirFirstCreateFlag = true; 

/// Auxiliary variable (initialization)
bool VertexDirFirstCreateFlag = true;

// Zero vector
vec2D ZeroVec2D = { 0, 0 };

// Unit vector
vec2D UnitVecR = { 1, 0 };

// Unit tensor
tensor2D UnitTensor = {1, 0, 0, 1};

// Zero tensor
tensor2D ZeroTensor = {0, 0, 0, 0};

/// Mesh structure
struct meshObject MeshObj;

// Material structure
/// 0th element - meaningless, 
/// 1st element - insulator type wall material, 
/// 2nd element - ideal conductor type wall material.
struct materialParameters MatWall[3];

/// Collection of materials
struct materialParameters MatParasList[MAX_MAT_TYPE];

// The number of materials
int MatParasListLen;


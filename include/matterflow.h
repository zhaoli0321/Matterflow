/*! \file  matterflow.h
 *
 *  \brief Main header file for the Matter flow project
 *
 *  \note  This header file contains general constants and data structures of
 *         Matterflow. It contains macros and data structure definitions; should not
 *         include function declarations here.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

//////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <malloc.h>


#ifndef __MATTERFLOW_HEADER__      /*-- allow multiple inclusions --*/
#define __MATTERFLOW_HEADER__      /**< indicate matterflow.h has been included before */


/*---------------------------*/
/*---  Macros definition  ---*/
/*---------------------------*/


#define MAX_VERTS 3000000 //500000;
#define MAX_TRKPT 10
#define MAX_NB_TRGS 20
#define MAX_MAT_TYPE 10


////////////////////////////////////////////////////// define bool, OFF, ON
#ifndef bool
#define bool int
#endif
#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif
#ifndef OFF
#define OFF 0
#endif
#ifndef ON
#define ON 1
#endif
////////////////////////////////////////////////////// define max, min
#ifndef max
#define max(a,b)  (((a)>(b)) ? (a):(b))
#endif
///////////
#ifndef min
#define min(a,b)  (((a)<(b)) ? (a):(b))
#endif
////////////////////////////////////////////////////// define PI
#ifndef PI
#define PI 3.1415926535897931
#endif

#define NullTriangle NULL


typedef struct vectorDimTwo
{
    double X;
    double Y;
} vec2D;


typedef struct tensorDimTwo
{
	double A[2][2];
} tensor2D;


typedef struct vertex
{
    /// Nodal index
    int Index;

    /// Marks if a node is "dead" (for future use)
    bool IsDead;

    /// Mark whether the node is on the up and down edges
    bool IsOnUpDownBoundary;
    bool IsOnUpBoundary;    // zhaoli
    bool IsOnDownBoundary;  // zhaoli

    /// Mark whether the node is on the left and right edges
    bool IsOnLeftRightBoundary;
    bool IsOnLeftBoundary;  // zhaoli
    bool IsOnRightBoundary; // zhaoli

    /// Nodal mass
    double Mass;

    /// Nodal velocity
    vec2D Velocity;

    /// Nodal positon
    vec2D Pos;

    ////////////////////////////////// Auxiliary variables
    /// Nodal force
    vec2D Force;

    /// Matter flow force (part of the nodal force)
    vec2D ForceForMatterFlow;

    /// displacement of one time step (Used to assist in calculating the internal energy increment.)
	vec2D DeltPos;
} *Vertex;


typedef struct triangle
{
    /// Triangle index
    int Index;

    /// Marks if a triangle cell is "dead" (for future use)
    bool IsDead;
    
    /// The three vertices of the triangle (conventionally arranged in a counterclockwise direction)
    Vertex Vertices[3];

    /// Three neighborhood triangles of the triangle
    struct triangle * NeighbourTrgs[3];

    /// Matter flow acceleration corresponding to each edge for the triangle
    double FlowAcc[3];

    /// Matter flow velocity corresponding to each edge for the triangle
    double FlowVelocity[3];

    /// Triangle mass
    double Mass;

    /// Material ID of the triangle
    int MaterialId;

    /// Density of the triangle
    double Density;

    /// Pressure of the triangle
	double Pressure;

    /// Sound velocity of the triangle
    double SoundVelocitySum;

    /// Viscosity coefficient
    double ViscCoeff;

    /// Scalar viscous stress
    double Viscos;  

    /// The sum of internal stresses in the fluid. 
    /// Including: pressure, elastic stress and viscous stress. But it does not include stress tensor in the magnetic field.
	tensor2D TotalFluidStress;

    /// Area of the triangle
    double Area;

    /// Heights of the triangle
    double Heights[3];

    /// Minimum height of the triangle
    double MinHeight;

    /// Three vertex coordinates (auxiliary variable, for periodic meshes)
	vec2D CycledPoses[3];

    /// Internal energy of the triangle
    double InternalEnergy;

    /// Temperature of the triangle
    double Temperature;
} *Triangle;


typedef enum boxBoundaryCondition
{
	BoxBoundaryCondition_Wall, // Solid Wall
	BoxBoundaryCondition_Cycle // Periodic
} BoxBoundaryCondition;


typedef struct meshObject
{
    /////////////////////////////////////////////////// Mesh
    /// The set of nodes
	struct vertex * Vertices;

    /// The number of nodes
    int VertsArrLen;

    /// The set of dead nodes
	int * DeadVertIdsArr;

    /// The number of dead nodes
    int DeadVertIdsArrLen;
    
    /// The set of triangles
	struct triangle * Trgs;

    /// The number of triangles
    int TrgsArrLen;

    /// The number of triangles
	int * DeadTrgIdsArr;

    /// The number of dead triangles
    int DeadTrgIdsArrLen;

    /////////////////////////////////////////////////// Boundary conditions
    BoxBoundaryCondition TopBottomBoundaryCondition;
	BoxBoundaryCondition LeftRightBoundaryCondition;
    
    /// Gravitational acceleration
    double GravityFactor;

    /// Calculate domain size
    double Width;
    double Height;
} *MeshObject;


typedef enum eOSTypeEnum { 
	EOSTypeEnum_Linear, 
	EOSTypeEnum_IdealGas
} EOSTypeEnum;


typedef struct materialParameters
{
    //////////////////////////////////////////////////////////////////////////////////////////// Material properties
    int MaterialNumber; // Corresponds to the material number in the user input file.
    //////////////////////////////////////////////////////////// Material name
    char MaterialName[101];
    //////////////////////////////////////////////////////////// EOS parameters
    EOSTypeEnum EOSType;     // material type
    /// Linear elastic fluid parameters
    double NormalDensity;    // room temperature/pressure density
    double Kai;              // bulk modulus
    double SpecificHeat;     // specific heat
    /// Polytropic  gas parameters
    double Gamma;            // polytropic index
    double MolMass;          // molar mass
    //////////////////////////////////////////////////////////// Viscosity parameter
    double MotionViscos;
} *Material;


typedef struct nodeFromFile
{
	/// Node coordinates
	vec2D Pos;
	/// Node velocity
	vec2D Vel;
	/// Adjacent triangles of the nodes
	int TrgIDs[MAX_NB_TRGS];
	int TrgIDsLen;
} *NodeFromFile;


typedef struct triangleFromFile
{
	/// Material ID
	int MaterialId;
    /// Pressure
    double Pressure;
    /// Temperature
    double Temperature;
	/// Three nodes of a triangle
	int NodeIDs[3];
	/// Three neighborhood triangles of the triangle
	int NeighbourTrgIDs[3];
} *TriangleFromFile;



/*
 *  declarations
 */

extern vec2D ZeroVec2D;
extern vec2D UnitVecR;
extern tensor2D UnitTensor;
extern tensor2D ZeroTensor;
extern struct meshObject MeshObj;
extern struct materialParameters MatParasList[MAX_MAT_TYPE];
extern struct materialParameters MatWall[3];
extern struct materialParameters MatParasList[MAX_MAT_TYPE];
extern int MatParasListLen;
extern double ROfGas; 

extern int HAVE_MF;     
extern int MF_Method;     
extern int TaylorGreen;     

extern char * modelInfoFileName;
extern double ViscCoeffTimes; 
extern double timeStepSafeFactor; 
extern int TimeNumbers;
extern double TimeEvolvedRecord;
extern int TimeCountForBitmap;
extern double TimeInterval;
extern double TimeEnd;
extern int TotBitmapOutputs;


extern bool CellDirFirstCreateFlag;
extern bool VertexDirFirstCreateFlag;
extern FILE * PhysQuantStatisticsFile;

extern int Iter;
extern double TotTime;
extern double MatterFlowTime;
extern clock_t Start_time;
extern clock_t End_time;




#endif /* end if for __MATTERFLOW_HEADER__ */

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

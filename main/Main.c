/*! \file  Main.c
 *
 *  \brief Main Function for Matter Flow Method
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"

// function statements
void PrintHelpInfo();
void SetParameter(int argc, char *argv[]);
void SetInitialOrDefaultValueOfGlobals();
void CreateOutPutFile();
void WriteOutputData();
void ImportModel(char *modelFile);
void Evolve();


// Main Function
void main(int argc, char *argv[])
{
	// Global variables initialization
	SetInitialOrDefaultValueOfGlobals();

	// Get command line parameters
	SetParameter(argc, argv);
		
	// Import a computational model
	ImportModel(modelInfoFileName);

	// Print basic information
	PrintInfo(); 
	PrintMeshScale();

	// Create output files
	CreateOutPutFile();

	// Evolution (numerical simulation)
	printf("Evolving: \r\n");
	Start_time = clock();
	Evolve();
	End_time = clock();
	TotTime = GetWallTime(Start_time, End_time);

	// Print Numerical Simulation Information
	PrintNSInformation();

	if(TaylorGreen){	
		TaylorGreenOutPutBottomEdgePhysics();
	// 	VolumeWeightedL1Norm();
	// 	TaylorGreenL2Norm();
	}	
}


void PrintHelpInfo()
{
	printf("\n");
	printf("-f : Model file.\n");
	printf("-mf: Is there the matter flow method (default %d), 1: yes, 0: no.\n", HAVE_MF);
	printf("-k : Times of viscosity coefficient (default %.4f).\n", ViscCoeffTimes);
	printf("-cs: C_safe is CFL constant (default %.4f).\n", timeStepSafeFactor);
	printf("-bn: How many pictures (or Tecplot files) to output (default %d).\n", TotBitmapOutputs);
	printf("\n");
}


// Read arguments from the command line
void SetParameter(int argc, char *argv[])
{
	/* Read arguments */
	int i = 1;
	while(i < argc)
	{	
		if(!strcmp(argv[i], "-f"))
		{
			modelInfoFileName = argv[i + 1];
		}
		if(!strcmp(argv[i], "-mf"))
		{
			HAVE_MF = atoi(argv[i + 1]);
		}	
		if(!strcmp(argv[i], "-bn"))
		{
			TotBitmapOutputs = atoi(argv[i + 1]);
		}	
		if(!strcmp(argv[i], "-k"))
		{
			ViscCoeffTimes = atof(argv[i + 1]);
		}	
		if(!strcmp(argv[i], "-cs"))
		{
			timeStepSafeFactor = atof(argv[i + 1]);
		}	
		// Print help information
		if(!strcmp(argv[i], "-help") || !strcmp(argv[i], "--help") || !strcmp(argv[i], "-h") || !strcmp(argv[i], "--h"))
		{
			printf("Options %s: \n", argv[i]);
			PrintHelpInfo();
			exit(-1);
		}		
		i += 2;
	}
}



void SetInitialOrDefaultValueOfGlobals()
{
	// Node array and triangle array initialization
	CreateVerticesArray(&MeshObj.Vertices);
	CreateTrgsArray(&MeshObj.Trgs);
	//////////////////////////////////////////// Physical Constants
	/// Gas universal constant
	ROfGas = 8.3149 * 1e-5; // 8.3149*J.mol^-1.k^-1
	//////////////////////////////////////////// Runtime control parameters initialization
	/// Records the time of system evolution for continued calculation.
	TimeEvolvedRecord = 0.0;
	/// The number of times the system has evolved. (it is equal to time divided by some interval)
	TimeNumbers = 0;
	/// Time interval used to control data output
	TimeInterval = 1.0;
	/// Used to control the end time of the program running
	TimeEnd = 1.0;
}

void CreateOutPutFile()
{
	PhysQuantStatisticsFile = MyOpenFile("Physical_Quantity_Statistics.txt", "w");
	{
		FILE * fp = PhysQuantStatisticsFile;
		fprintf(fp, "The statistics of physical quantities are shown in the table below.");
		fprintf(fp, "\n\n");
		fprintf(fp, "     Time(us)   \t");
		/// List the physical quantities that you want to output in turn
		{
			fprintf(fp, "   Mass(g)    ");
			fprintf(fp, "   Momentum_x(g)");
			fprintf(fp, "   Momentum_y(g)");
			fprintf(fp, "   Energy(100KJ)");
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
}

// Write statistical physical quantities
void WriteOutputData()
{
	int i;
	// Count and output various physical quantity statistics that you want to output in turn
	double totMass = 0.0, totMomentumX = 0.0, totMomentumY = 0.0;
	double totEnergy = 0.0;
	for (i = 0; i < MeshObj.TrgsArrLen; i++)
	{
		Triangle trg = &MeshObj.Trgs[i];
		///////////////////////////////////////
		totMass += trg->Mass;
	}
	for (i = 0; i < MeshObj.VertsArrLen; i++)
	{
		Vertex vert = &MeshObj.Vertices[i];
		//////////////////////////////
		totMomentumX += vert->Mass * vert->Velocity.X;
	}
	for (i = 0; i < MeshObj.VertsArrLen; i++)
	{
		Vertex vert = &MeshObj.Vertices[i];
		//////////////////////////////
		totMomentumY += vert->Mass * vert->Velocity.Y;
	}
	for (i = 0; i < MeshObj.TrgsArrLen; i++)
	{
		Triangle trg = &MeshObj.Trgs[i];
		//////////////////////////////
		if (trg->MaterialId > 0)
		{
			totEnergy += trg->InternalEnergy;
		}
	}
	for (i = 0; i < MeshObj.VertsArrLen; i++)
	{
		Vertex vert = &MeshObj.Vertices[i];
		//////////////////////////////
		totEnergy += 0.5 * vert->Mass * CalcLengthSquar(vert->Velocity);
		totEnergy += vert->Mass * vert->Pos.Y * MeshObj.GravityFactor; // add gravitational potential energy
	}
	/// Output one by one
	fprintf(PhysQuantStatisticsFile, "%15e", TimeEvolvedRecord);
	fprintf(PhysQuantStatisticsFile, "\t%15e", totMass);
	fprintf(PhysQuantStatisticsFile, "\t%15e", totMomentumX);
	fprintf(PhysQuantStatisticsFile, "\t%15e", totMomentumY);
	fprintf(PhysQuantStatisticsFile, "\t%15e", totEnergy);
	/////////////////////////////////////////////////////////
	fprintf(PhysQuantStatisticsFile, "\n");
	fflush(PhysQuantStatisticsFile);
}

// Import physical conditions and meshes from modeling files
void ImportModel(char * modelFile)
{
	printf("Loading compute model information from file \"%s\" \r\n\r\n", modelFile);
	MeshObjectFromLoadModel(modelFile);
	printf("Compute model information loaded. \r\n\r\n");
}


void Evolve()
{
	clock_t t1, t2, temp_time;
	double deltT, Time;

	/// Set dependent variables
	SetAllDependentVariablesOfTrgs();

	//////////////////////////////////////////////////////////////////////////////////////// Output initial state data file and Tecplot
	printf("%3d: Iter = %6d, t = %-.6f\n", TimeCountForBitmap, Iter, TimeEvolvedRecord);
	WriteOutputData();
	WriteTecPlot();

	////////////////////////////////////////// Data preparation ///////////////////////////////////////////
	/// Set velocity, position boundary conditions
	VelocityPosEdgeCondition();
	/// Set gravity 
	SetGravity();
	/// Calculate the force of the node 
	CalculateVertexForce();
	/// Set force boundary condition
	ForceEdgeCondition();
	/// Determine the time-step  
	deltT = DetermineDeltT();

	////////////////////////////////////////// T = T + \Delta T ///////////////////////////////////////////
	while (TimeEvolvedRecord <= TimeEnd)
	{
		Iter++;

	    /// Update the physical quantities at the next moment
		DynamicEvolve(deltT);
	    /// Update dependent variables
	    SetAllDependentVariablesOfTrgs();
	    /// Matter flow evolution (alleviate the checkerboard oscillations of physical quantities)
		if(HAVE_MF){
			t1 = clock();
			MatterFlowEvolve(deltT);
			t2 = clock();
			MatterFlowTime += GetWallTime(t1, t2);
		}

		//////////////////////////////////////////// Prepare for the next time-step
	    /// Set velocity, position boundary conditions
		VelocityPosEdgeCondition();
		/// Set gravity 
		SetGravity();
		/// Calculate the force of the node 
	    CalculateVertexForce();
		/// Set force boundary condition
	    ForceEdgeCondition();
	    /// Determine the time-step  
		deltT = DetermineDeltT();

	    //////////////////////////////////////////// Post-processing
	    TimeEvolvedRecord += deltT;
		printf("\r  Iter = %6d, t = %10.6e, deltT = %10.6e", Iter, TimeEvolvedRecord, deltT);
	    /////////////////////////////////////////// output files
	    if (TimeEvolvedRecord / (TimeEnd / TotBitmapOutputs) > TimeCountForBitmap + 1)
	    {
	        TimeCountForBitmap++;
	        //////////////////////////////// output
			temp_time = clock();
			Time = GetWallTime(Start_time, temp_time);
			printf("\r%3d: Iter = %6d, t = %10.6e, deltT = %10.6e, TimeUsed=%.2fs\n", TimeCountForBitmap, Iter, TimeEvolvedRecord, deltT, Time);
			WriteTecPlot();
	    }
	    if (TimeEvolvedRecord / TimeInterval > TimeNumbers + 1)
	    {
	        TimeNumbers++;
	        //////////////////////////////////////////// Write physical quantities to data file
	        WriteOutputData();
	    }
	}
	
}








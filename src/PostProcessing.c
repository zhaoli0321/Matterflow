/*! \file  PostProcessing.c
 *
 *  \brief Define some PostProcessing functions
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"


/* Create directory: determine whether the directory exists, if not, create the directory */
void CreateDirectory(char *dirName)
{			
	int isExistDir = access(dirName, 0); // 0(exist), -1(does not exist)
	if(isExistDir == -1){
		char mkdir[128];
		sprintf(mkdir, "mkdir %s", dirName);
		//printf("dirName = %s\n", dirName);
		system(mkdir); // create directory
	}
}


// Numerical Simulation Information
void PrintNSInformation()
{
    printf("\r\n\r\n");
    printf("                            Numerical Simulation Information                                \r\n");
    printf("--------------------------------------------------------------------------------------------\r\n");
    printf("The number of triangles: %d, the number of vertices: %d, Iter = %d\r\n", MeshObj.TrgsArrLen, MeshObj.VertsArrLen, Iter);
    printf("Simulation Time: \t\t %.4f\r\n", TotTime);
    printf("MatterFlow Time: \t\t %.4f (%.2f%%)\r\n", MatterFlowTime, 100 * MatterFlowTime / TotTime);
    printf("--------------------------------------------------------------------------------------------\r\n\r\n");
}


void PrintInfo()
{
	printf("HAVE_MF = %d, TotBitmapOutputs = %d \n\n", HAVE_MF, TotBitmapOutputs);
	printf("Times of viscosity coefficient %.4f, C_safe = %.4f\n\n", ViscCoeffTimes, timeStepSafeFactor);
	printf("Domain Omega = [0, %.2f] x [0, %.2f]\n", MeshObj.Width, MeshObj.Height); // print domain size
	printf("\n");
}


//// Output mesh size
void PrintMeshScale()
{
	int i, k;
	double side, minSideLen = 1e10, maxSideLen = 0;
	Triangle trg;
	///////////////////////////////////////
	for (k = 0; k < MeshObj.TrgsArrLen; k++)
	{
		trg = &MeshObj.Trgs[k];
		///////
		for (i = 0; i < 3; i++)
		{
			side = CalcLength(Vec2DSub(trg->CycledPoses[(i + 1) % 3], trg->CycledPoses[i]));
			minSideLen = min(side, minSideLen);
			maxSideLen = max(side, maxSideLen);
		}
	}
	printf("MinMeshScale = %f, MaxMeshScale = %f, AvgMeshScale = %f \r\n\r\n", minSideLen, maxSideLen, (minSideLen + maxSideLen) / 2.0);				
}


// Write cell data and node data in TecPlot format (TecPlot software)
// The output directory is TecPlot/cPhysicsFile (cell data)
// The output directory is TecPlot/vPhysicsFile (node data)
void WriteTecPlot()
{
    /////////////////////// declare variables ////////////////////////////
    char        fileName[256];
	char        rootDir[64] = "TecPlot";
	char        subDir[64] = "cPhysicsFile";
	char        subsubDir[128];
	char 		currentDir[512];
	char        prefixFileName[64] = "cPhysics";
    int         *elem = NULL, *vMap = NULL;
    double      *node = NULL;
    int         v, c, vSize = MeshObj.VertsArrLen, cSize = MeshObj.TrgsArrLen;
    ////////////////////////////////////////////////////////////

    /*------------------------------------------------------------------------------------------------*/
    // cell data
    /*------------------------------------------------------------------------------------------------*/
    /* Determine whether the directory exists, if not, create the directory */
    if(CellDirFirstCreateFlag){
        CreateDirectory(rootDir);
        sprintf(currentDir, "%s/%s", rootDir, subDir);
        CreateDirectory(currentDir);
        // subdirectory
        {
            sprintf(subsubDir, "%s", currentDir);
            int isExistSubDir = access(subsubDir, 0);
            if(isExistSubDir == -1){
                char mkdir[512];
                sprintf(mkdir, "mkdir %s", subsubDir);
                system(mkdir);
            }else{
                char mkdir[512], rm[512];
                sprintf(mkdir, "mkdir %s", subsubDir);
                sprintf(rm, "rm -r %s", subsubDir);

                system(rm); // first delete the directory
                system(mkdir); // then create the directory
            }
        }
        CellDirFirstCreateFlag = false;								
    }
    /////////////
    sprintf(currentDir, "%s/%s", rootDir, subDir);
    sprintf(fileName, "%s/%s_%d_%d.dat", currentDir, prefixFileName, TimeCountForBitmap, Iter);


// Count the number of vacuum elements
int VacuumNEM = 0;
#if Omit_Vacuum_Element     
    for(c = 0; c < cSize; c++)
    {
        if(MeshObj.Trgs[c].MaterialId <= 0) VacuumNEM++;
    }
#endif 


    /* write tecPlot */
    FILE  *fp = NULL;
    // tecplot
    fp = fopen(fileName, "w");

    // write head file in TecPlot format 
    fprintf(fp, "TITLE = \"2D: Cells Data, Time = %e s\"\n", TimeEvolvedRecord);
    fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Pressure\", \"Density\", \"Temperature\", \"MaterialID\", \"Mass\", \"InternalEnergy\"\n");
    fprintf(fp, "ZONE T = \"%e\", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, VARLOCATION=([3-8]=CELLCENTERED), ZONETYPE=FETRIANGLE\n",
                                TimeEvolvedRecord, vSize, cSize - VacuumNEM);
    

    // write node coordinates and variables
    // 1) x_coord
    double X, Y;
    for  (v = 0; v < vSize; v++)
    {
        X = MeshObj.Vertices[v].Pos.X;
        fprintf(fp, "%f\n", X);
    }

    // 2) y_coord
    for  (v = 0; v < vSize; v++)
    {
        Y = MeshObj.Vertices[v].Pos.Y;
        fprintf(fp, "%f\n", Y);
    }

    // 3) Pressure
    for(c = 0; c < cSize; c++)
    {
#if Omit_Vacuum_Element        
        if(MeshObj.Trgs[c].MaterialId <= 0) continue;
#endif          
        fprintf(fp, "%.20f\n", MeshObj.Trgs[c].Pressure);
    }

    // 4) Density
    for(c = 0; c < cSize; c++)
    {
#if Omit_Vacuum_Element        
        if(MeshObj.Trgs[c].MaterialId <= 0) continue;
#endif          
        fprintf(fp, "%.20f\n", MeshObj.Trgs[c].Density);
    }

    // 5) Temperature
    for(c = 0; c < cSize; c++)
    {
#if Omit_Vacuum_Element        
        if(MeshObj.Trgs[c].MaterialId <= 0) continue;
#endif          
        fprintf(fp, "%.20f\n", MeshObj.Trgs[c].Temperature);
    }

    // 6) MatrialID
    for(c = 0; c < cSize; c++)
    {
#if Omit_Vacuum_Element        
        if(MeshObj.Trgs[c].MaterialId <= 0) continue;
#endif          
        fprintf(fp, "%d\n", MeshObj.Trgs[c].MaterialId);
    }

    // 7) Mass
    for(c = 0; c < cSize; c++)
    {
#if Omit_Vacuum_Element        
        if(MeshObj.Trgs[c].MaterialId <= 0) continue;
#endif          
        fprintf(fp, "%.20f\n", MeshObj.Trgs[c].Mass);
    }
    
    // 8) InternalEnergy
    for(c = 0; c < cSize; c++)
    {
#if Omit_Vacuum_Element        
        if(MeshObj.Trgs[c].MaterialId <= 0) continue;
#endif      
        fprintf(fp, "%.20f\n", MeshObj.Trgs[c].InternalEnergy);
    }
    
    // write elem
    for(c = 0; c < cSize; c++)
    {
#if Omit_Vacuum_Element        
        if(MeshObj.Trgs[c].MaterialId <= 0) continue;
#endif 
        fprintf(fp, "%d %d %d\n", MeshObj.Trgs[c].Vertices[0]->Index + 1, \
                    MeshObj.Trgs[c].Vertices[1]->Index + 1, MeshObj.Trgs[c].Vertices[2]->Index + 1);
    }
    ////
    fclose(fp);



    /*------------------------------------------------------------------------------------------------*/
    // Node data
    /*------------------------------------------------------------------------------------------------*/
    /////////////////////// declare variables ////////////////////////////
	char			vfileName[256];
	char        	vrootDir[64] = "TecPlot";
	char        	vsubDir[64] = "vPhysicsFile";
	char        	vsubsubDir[128];
	char 			vcurrentDir[512];
	char			vprefixFileName[30] = "vPhysics";
    Vertex          vert;

    /* Determine whether the directory exists, if not, create the directory */
    if(VertexDirFirstCreateFlag){
        CreateDirectory(vrootDir);
        sprintf(vcurrentDir, "%s/%s", vrootDir, vsubDir);
        CreateDirectory(vcurrentDir);
        // subdirectory
        {
            sprintf(vsubsubDir, "%s", vcurrentDir);
            int isExistSubDir = access(vsubsubDir, 0);
            if(isExistSubDir == -1){
                char mkdir[512];
                sprintf(mkdir, "mkdir %s", vsubsubDir);
                system(mkdir);
            }else{
                char mkdir[512], rm[512];
                sprintf(mkdir, "mkdir %s", vsubsubDir);
                sprintf(rm, "rm -r %s", vsubsubDir);

                system(rm);     // first delete the directory
                system(mkdir);  // then create the directory
            }
        }
        VertexDirFirstCreateFlag = false;								
    }
    /////////////
    sprintf(vcurrentDir, "%s/%s", vrootDir, vsubDir);
    sprintf(vfileName, "%s/%s_%d_%d.dat", vcurrentDir, vprefixFileName, TimeCountForBitmap, Iter);
    
    
    /* write tecPlot */
    // tecplot
    fp = fopen(vfileName, "w");
    
    // write head file in TecPlot format 
    fprintf(fp, "TITLE = \"2D: Vertices Data, Time = %e s\"\n", TimeEvolvedRecord);
    fprintf(fp, "VARIABLES = \"X\", \"Y\", "
                "\"Mass\", " 
                "\"Velocity_X\", "
                "\"Velocity_Y\", "
                "\"Velocity\", "
                "\"Force_X\", "
                "\"Force_Y\", "
                "\"ForceForMatterFlow_X\", "
                "\"ForceForMatterFlow_Y\"\n");
    fprintf(fp, "ZONE T = \"%e\", NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n", TimeEvolvedRecord, vSize, cSize - VacuumNEM);
    
    ////////////////////////////////////////////////////////////
    // write node coordinates and variables
    for  (v = 0; v < vSize; v++)
    {
        vert = &MeshObj.Vertices[v];
        X = vert->Pos.X;
        Y = vert->Pos.Y;

        fprintf(fp, "%f %f "
                    "%.20f "
                    "%.20f "
                    "%.20f "
                    "%.20f "
                    "%.20f "
                    "%.20f "
                    "%.20f "
                    "%.20f \n",

                    X, Y,                               // x, y
                    vert->Mass, 		                // Mass
                    vert->Velocity.X, 	                // Velocity_X
                    vert->Velocity.Y,	                // Velocity_Y
                    sqrt(vert->Velocity.X*vert->Velocity.X+vert->Velocity.Y*vert->Velocity.Y), // Velocity size
                    vert->Force.X,		                // Force_X
                    vert->Force.Y, 		                // Force_Y
                    vert->ForceForMatterFlow.X, 		// ForceForMatterFlow_X
                    vert->ForceForMatterFlow.Y); 		// ForceForMatterFlow_Y
    }
    
    // write elem
    for(c = 0; c < cSize; c++)
    {
#if Omit_Vacuum_Element        
        if(MeshObj.Trgs[c].MaterialId <= 0) continue;
#endif        
        fprintf(fp, "%d %d %d\n", MeshObj.Trgs[c].Vertices[0]->Index + 1, \
                    MeshObj.Trgs[c].Vertices[1]->Index + 1, MeshObj.Trgs[c].Vertices[2]->Index + 1);
    }
    ////
    fclose(fp);
}



/*------------------------------------------------------------------------------------------------*/
// Computational error for Taylor-Green problem, 2021.11.08
/*------------------------------------------------------------------------------------------------*/
// Pressure analytical solution
double PressureAnalyticalSolution_TaylorGreen(double xc, double yc)
{
    double rho = 1;
    return (10 + rho / 4.0 *( cos(2*PI*xc) + cos(2*PI*yc) ) );
}


// Velocity Analytical Solution
void VelocityAnalyticalSolution_TaylorGreen(double x, double y, double *vx, double *vy)
{
    *vx =   sin(PI*x)*cos(PI*y);
    *vy = - cos(PI*x)*sin(PI*y);
}


// For Taylor-Green problem, the error of calculation velocity and pressure
void TaylorGreenL2Norm()
{
    int         v = 0, c = 0, vSize = MeshObj.VertsArrLen, cSize = MeshObj.TrgsArrLen;
    double      x = 0, y = 0, xc = 0, yc = 0;
    double      VelL2 = 0, PreL2 = 0;
    double      vx_exact, vy_exact, vx, vy, p_exact, p;

    /////// For node 
    for  (v = 0; v < vSize; v++)
    {
        x = MeshObj.Vertices[v].Pos.X;
        y = MeshObj.Vertices[v].Pos.Y;
        
        // Computational velocity analytical solution
        VelocityAnalyticalSolution_TaylorGreen(x, y, &vx_exact, &vy_exact);
        
        vx = MeshObj.Vertices[v].Velocity.X;
        vy = MeshObj.Vertices[v].Velocity.Y;

        VelL2 += (vx - vx_exact)*(vx - vx_exact) + (vy - vy_exact)*(vy - vy_exact);
    }
    VelL2 /= sqrt(vSize);


    /////// For cell
    for(c = 0; c < cSize; c++)
    {
        xc = (MeshObj.Trgs[c].Vertices[0]->Pos.X + MeshObj.Trgs[c].Vertices[1]->Pos.X + MeshObj.Trgs[c].Vertices[2]->Pos.X) / 3.0;
        yc = (MeshObj.Trgs[c].Vertices[0]->Pos.Y + MeshObj.Trgs[c].Vertices[1]->Pos.Y + MeshObj.Trgs[c].Vertices[2]->Pos.Y) / 3.0;
        
        p_exact = PressureAnalyticalSolution_TaylorGreen(xc, yc);
        p = MeshObj.Trgs[c].Pressure;

        PreL2 += (p_exact - p)*(p_exact - p);
    }
    PreL2 /= sqrt(cSize);

    printf("Print TaylorGreenL2Norm\n");
    printf("VelL2 = %e\n", VelL2);
    printf("PreL2 = %e\n", PreL2);
}


// Volume-Weighted L1 error of calculation speed and pressure for Taylor-Green problem
void VolumeWeightedL1Norm()
{
    int         v = 0, c = 0, vSize = MeshObj.VertsArrLen, cSize = MeshObj.TrgsArrLen;
    double      x = 0, y = 0, xc = 0, yc = 0;
    double      VelL1 = 0, PreL1 = 0, TotalNodalVolume = 0, TotalElementVolume = 0;
    double      vx_exact, vy_exact, vx, vy, p_exact, p, elem_vol, *node_vol = NULL;

    node_vol = (double *) calloc(vSize, sizeof(double));
    /////// Calculate the volume of the node control volum
    for(c = 0; c < cSize; c++)
    {
        xc = (MeshObj.Trgs[c].Vertices[0]->Pos.X + MeshObj.Trgs[c].Vertices[1]->Pos.X + MeshObj.Trgs[c].Vertices[2]->Pos.X) / 3.0;
        yc = (MeshObj.Trgs[c].Vertices[0]->Pos.Y + MeshObj.Trgs[c].Vertices[1]->Pos.Y + MeshObj.Trgs[c].Vertices[2]->Pos.Y) / 3.0;
        
        // printf("trg[%d], {(%.2f, %0.2f), (%.2f, %0.2f), (%.2f, %0.2f)}\n", c, MeshObj.Trgs[c].Vertices[0]->Pos.X, MeshObj.Trgs[c].Vertices[0]->Pos.Y, 
        // MeshObj.Trgs[c].Vertices[1]->Pos.X, MeshObj.Trgs[c].Vertices[1]->Pos.Y, MeshObj.Trgs[c].Vertices[2]->Pos.X, MeshObj.Trgs[c].Vertices[2]->Pos.Y);

        p_exact = PressureAnalyticalSolution_TaylorGreen(xc, yc);
        p = MeshObj.Trgs[c].Pressure;
        elem_vol = MeshObj.Trgs[c].Area;
        
        // Calculate the volume of the node control volume
        node_vol[MeshObj.Trgs[c].Vertices[0]->Index] += elem_vol / 3.0;
        node_vol[MeshObj.Trgs[c].Vertices[1]->Index] += elem_vol / 3.0;
        node_vol[MeshObj.Trgs[c].Vertices[2]->Index] += elem_vol / 3.0;

        if (MeshObj.Trgs[c].MaterialId > 0){
            // Calculate Volume-Weighted L1 Error of Pressure
            TotalElementVolume += elem_vol;
            PreL1 += fabs(p_exact - p)*elem_vol;
        }
    }
    PreL1 /= TotalElementVolume;

    for  (v = 0; v < vSize; v++)
    {
        x = MeshObj.Vertices[v].Pos.X;
        y = MeshObj.Vertices[v].Pos.Y;
        
        // Computational velocity analytical solution
        VelocityAnalyticalSolution_TaylorGreen(x, y, &vx_exact, &vy_exact);
        
        vx = MeshObj.Vertices[v].Velocity.X;
        vy = MeshObj.Vertices[v].Velocity.Y;

        TotalNodalVolume += node_vol[v];
        VelL1 += sqrt( (vx - vx_exact)*(vx - vx_exact) + (vy - vy_exact)*(vy - vy_exact) )*node_vol[v];
    }
    VelL1 /= TotalNodalVolume;

    // free
    free(node_vol); node_vol = NULL;


    // output
    printf("TaylorGreen Volume-Weighted L1 error Norm:\n");
    printf("VelL1 = %e, TotalNodalVolume   = %e\n", VelL1, TotalNodalVolume);
    printf("PreL1 = %e, TotalElementVolume = %e\n", PreL1, TotalElementVolume);
}


void TaylorGreenOutPutBottomEdgePhysics()
{
    int         v = 0, c = 0, vSize = MeshObj.VertsArrLen, cSize = MeshObj.TrgsArrLen;
    double      x = 0, y = 0, xc = 0, yc = 0;
    double      vx_exact, vy_exact, vx, vy, p_exact, p;
    FILE        *fp = NULL;
    char        vfileName[128] = "";
    char        cfileName[128] = "";

    /////// For node
    sprintf(vfileName, "%s.txt", "Velocity_x_BottomEdge");
    fp = fopen(vfileName, "w");

    fprintf(fp, "x\tvelocity_x\texact\n");
    for  (v = 0; v < vSize; v++)
    {
        if (!MeshObj.Vertices[v].IsOnDownBoundary) continue;
        
        x = MeshObj.Vertices[v].Pos.X;
        y = MeshObj.Vertices[v].Pos.Y;
        
        // Computational velocity analytical solution
        VelocityAnalyticalSolution_TaylorGreen(x, y, &vx_exact, &vy_exact);
        
        vx = MeshObj.Vertices[v].Velocity.X;
        vy = MeshObj.Vertices[v].Velocity.Y;

        fprintf(fp, "%e\t%e\t%e\n", x, vx, vx_exact);
    }
    ////
    fclose(fp);

       
    /////// For cell
    // Generates an array of triangle numbers with adjacent vertices
    int *vnbtrgmark     = (int *) calloc(vSize, sizeof(int));
    int (*vnbtrg)[25]   = (int(*)[25]) calloc(vSize*25, sizeof(int));
    // int (*vnbtrg)[25]   = (int *) calloc(vSize*25, sizeof(int));
    int k, t;

	/// Traversal triangular element
	for (c = 0; c < cSize; c++)
	{
		for (v = 0; v < 3; v++ )
		{
			k = MeshObj.Trgs[c].Vertices[v]->Index; // nodal index
			t = vnbtrgmark[k];
			vnbtrg[k][t] = c; // Store the index of adjacent triangle elements
			vnbtrgmark[k]++;  // Count the number of adjacent triangle elements		
		}
	}

    sprintf(cfileName, "%s.txt", "Pressure_BottomEdge");
    fp = fopen(cfileName, "w");

    fprintf(fp, "x\tpressure\texact\n");
    for (v = 0; v < vSize; v++)
    {
        if (!MeshObj.Vertices[v].IsOnDownBoundary) continue;
        
        x = MeshObj.Vertices[v].Pos.X;
        y = MeshObj.Vertices[v].Pos.Y;
        
        // Calculate the pressure analytical solution
        p_exact = PressureAnalyticalSolution_TaylorGreen(x, y);
        
        p = 0.0;
        for(c = 0; c < vnbtrgmark[v]; c++){
            p += MeshObj.Trgs[ vnbtrg[v][c] ].Pressure;
        }
        p /= vnbtrgmark[v];
        // printf("%f %f %f\n", x, p, p_exact);

        fprintf(fp, "%e\t%e\t%e\n", x, p, p_exact);
    }
    ////
    fclose(fp);

    // 
    free(vnbtrgmark);
    free(vnbtrg);
}
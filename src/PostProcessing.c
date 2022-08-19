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

// internal function declarations
int system(char *);
int access(char *, int );

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
	printf("Times of viscosity coefficient %.4f\n\n", ViscCoeffTimes);
	printf("Times of heat diffusion coefficient %.4f\n\n", TDCoeffTimes);
	printf("Time-step factor C_safe = %.4f\n\n", timeStepSafeFactor);
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
    double      VelL2 = 0, PreL2 = 0, DenL2 = 0;
    double      vx_exact, vy_exact, vx, vy, p_exact, p, rho, rho_exact = 1;

    /// For node 
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
    VelL2 = sqrt(VelL2/vSize);


    /// For cell
    for(c = 0; c < cSize; c++)
    {
        xc = (MeshObj.Trgs[c].Vertices[0]->Pos.X + MeshObj.Trgs[c].Vertices[1]->Pos.X + MeshObj.Trgs[c].Vertices[2]->Pos.X) / 3.0;
        yc = (MeshObj.Trgs[c].Vertices[0]->Pos.Y + MeshObj.Trgs[c].Vertices[1]->Pos.Y + MeshObj.Trgs[c].Vertices[2]->Pos.Y) / 3.0;
        
        // Computational pressure analytical solution
        p_exact = PressureAnalyticalSolution_TaylorGreen(xc, yc);
        p = MeshObj.Trgs[c].Pressure;
        rho = MeshObj.Trgs[c].Density;

        PreL2 += (p_exact - p)*(p_exact - p);
        DenL2 += (rho_exact - rho)*(rho_exact - rho);
    }
    PreL2 = sqrt(PreL2/cSize);
    DenL2 = sqrt(DenL2/cSize);

    printf("Print TaylorGreenL2Norm\n");
    printf("Density  L2 = %e\n", DenL2);
    printf("Pressure L2 = %e\n", PreL2);
    printf("Velocity L2 = %e\n", VelL2);
}




/*------------------------------------------------------------------------------------------------*/
// Computational error for Gresho vortex problem, 2022.08.19
/*------------------------------------------------------------------------------------------------*/
void VelocityExactForGresho(double x, double y, double *vx, double *vy)
{
   // Vortex center is (xc, yc)
   double xc = MeshObj.Width / 2; 
   double yc = MeshObj.Width / 2; 
   x = x - xc;
   y = y - yc;
   const double r = sqrt(x * x + y * y);

   if (r < 0.2)
   {
      *vx =  5.0 * y;
      *vy = -5.0 * x;
   }
   else if (r < 0.4)
   {
      *vx =  2.0 * y / r - 5.0 * y;
      *vy = -2.0 * x / r + 5.0 * x;
   }
   else { 
      *vx = 0.0; 
      *vy = 0.0; 
   }
}

double PressureExactForGresho(double x, double y)
{
    // Vortex center is (xc, yc)
    double xc = MeshObj.Width / 2; 
    double yc = MeshObj.Width / 2; 
    const double rsq = (x-xc) * (x-xc) + (y-yc) * (y-yc), r = sqrt(rsq);
    if (r < 0.2)
    {
        return (5.0 + 25.0 / 2.0 * rsq);
    }
    else if (r < 0.4)
    {
        const double t1 = 9.0 - 4.0 * log(0.2) + 25.0 / 2.0 * rsq;
        const double t2 = 20.0 * r - 4.0 * log(r);
        return (t1 - t2);
    }
    else { return (3.0 + 4.0 * log(2.0)); }
}


/*
* L2 norm for calculating velocity and pressure
*/
void ComputeL2normForGresho()
{
    int         v = 0, c = 0, vSize = MeshObj.VertsArrLen, cSize = MeshObj.TrgsArrLen;
    double      x = 0, y = 0, xc = 0, yc = 0;
    double      VelL2 = 0, PreL2 = 0, DenL2 = 0;
    double      vx_exact, vy_exact, vx, vy, p_exact, p, rho_exact = 1, rho;

    /// Node variables
    for  (v = 0; v < vSize; v++)
    {
        x = MeshObj.Vertices[v].Pos.X;
        y = MeshObj.Vertices[v].Pos.Y;

        // Analytical solution of computational velocity
        VelocityExactForGresho(x, y, &vx_exact, &vy_exact);
        
        vx = MeshObj.Vertices[v].Velocity.X;
        vy = MeshObj.Vertices[v].Velocity.Y;

        VelL2 += (vx - vx_exact)*(vx - vx_exact) + (vy - vy_exact)*(vy - vy_exact);
    }

    VelL2 = sqrt(VelL2/vSize);

    /// Element variables
    for(c = 0; c < cSize; c++)
    {
        xc = (MeshObj.Trgs[c].Vertices[0]->Pos.X + MeshObj.Trgs[c].Vertices[1]->Pos.X + MeshObj.Trgs[c].Vertices[2]->Pos.X) / 3.0;
        yc = (MeshObj.Trgs[c].Vertices[0]->Pos.Y + MeshObj.Trgs[c].Vertices[1]->Pos.Y + MeshObj.Trgs[c].Vertices[2]->Pos.Y) / 3.0;
        
        // Analytical solution of computational pressure
        p_exact = PressureExactForGresho(xc, yc);
        p = MeshObj.Trgs[c].Pressure;

        PreL2 += (p_exact - p)*(p_exact - p);

        // density
        rho = MeshObj.Trgs[c].Density;
        DenL2 += (rho_exact - rho)*(rho_exact - rho);
    }

    PreL2 = sqrt(PreL2/cSize);
    DenL2 = sqrt(DenL2/cSize);

    // print information
    printf("Print Gresho L2 Norm:\n");
    printf("Density  L2 = %e\n", DenL2);
    printf("Pressure L2 = %e\n", PreL2);
    printf("Velocity L2 = %e\n", VelL2);
}

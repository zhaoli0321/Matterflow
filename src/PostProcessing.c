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

/*------------------------------------------------------------------------------------------------*/
// 后处理文件，打印模拟时间
/*------------------------------------------------------------------------------------------------*/

/* 创建目录: 判断目录是否存在，若不存在，则创建目录 */
void CreateDirectory(char *dirName)
{			
	int isExistDir = access(dirName, 0);// 0(存在), -1(不存在)
	if(isExistDir == -1){
		char mkdir[128];
		sprintf(mkdir, "mkdir %s", dirName);
		//printf("dirName = %s\n", dirName);
		system(mkdir);
	}
}



void PrintTime()
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
	printf("HAVE_MF = %d, MF_Method = %d, TotBitmapOutputs = %d \n\n", HAVE_MF, MF_Method, TotBitmapOutputs);
	printf("Times of viscosity coefficient %.4f, C_safe = %.4f\n\n", ViscCoeffTimes, timeStepSafeFactor);
	printf("Domain Omega = [0, %.2f] x [0, %.2f]\n", MeshObj.Width, MeshObj.Height); // 打印区域大小
	printf("\n");
}


//// 输出网格尺寸
void PrintMeshScale()
{
	int i, k;
	double side, minSideLen = 1e10, maxSideLen = 0;
	Triangle trg;
	///////////////////////////////////////
	for (k = 0; k < MeshObj.TrgsArrLen; k++)
	{
		trg = &MeshObj.Trgs[k];
		if (trg->IsDead) continue;
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


// 写TecPlot格式的单元数据和节点数据 
void WriteTecPlotCellDataFromNormalGrid()
{
    /////////////////////// 声明变量 ////////////////////////////
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
    // 单元数据
    /*------------------------------------------------------------------------------------------------*/
    /* 判断目录是否存在，若不存在，则创建目录 */
    if(CellDirFirstCreateFlag){
        CreateDirectory(rootDir);
        sprintf(currentDir, "%s/%s", rootDir, subDir);
        CreateDirectory(currentDir);
        // 子目录
        {
            sprintf(subsubDir, "%s", currentDir);
            int isExistSubDir = access(subsubDir, 0);// 0(存在), -1(不存在)
            if(isExistSubDir == -1){
                char mkdir[512];
                sprintf(mkdir, "mkdir %s", subsubDir);
                system(mkdir);
            }else{
                char mkdir[512], rm[512];
                sprintf(mkdir, "mkdir %s", subsubDir);
                sprintf(rm, "rm -r %s", subsubDir);

                system(rm);//先delete
                system(mkdir);// 再mkdir
            }
        }
        CellDirFirstCreateFlag = false;								
    }
    /////////////
    sprintf(currentDir, "%s/%s", rootDir, subDir);
    sprintf(fileName, "%s/%s_%d_%d.dat", currentDir, prefixFileName, TimeCountForBitmap, Iter);


// 统计真空单元数
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

    // write tecplot head file
    fprintf(fp, "TITLE = \"2D: Cells Data, Time = %e s\"\n", TimeEvolvedRecord);
    fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Pressure\", \"Density\", \"Temperature\", \"MaterialID\", \"Mass\", \"InternalEnergy\"\n");
    fprintf(fp, "ZONE T = \"%e\", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, VARLOCATION=([3-8]=CELLCENTERED), ZONETYPE=FETRIANGLE\n",
                                TimeEvolvedRecord, vSize, cSize - VacuumNEM);
    

    // write node and variable
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
    // 节点数据
    /*------------------------------------------------------------------------------------------------*/
    /////////////////////// 声明变量 ////////////////////////////
	char			vfileName[256];
	char        	vrootDir[64] = "TecPlot";
	char        	vsubDir[64] = "vPhysicsFile";
	char        	vsubsubDir[128];
	char 			vcurrentDir[512];
	char			vprefixFileName[30] = "vPhysics";
    Vertex          vert;

    /* 判断目录是否存在，若不存在，则创建目录 */
    if(VertexDirFirstCreateFlag){
        CreateDirectory(vrootDir);
        sprintf(vcurrentDir, "%s/%s", vrootDir, vsubDir);
        CreateDirectory(vcurrentDir);
        // 子目录
        {
            sprintf(vsubsubDir, "%s", vcurrentDir);
            int isExistSubDir = access(vsubsubDir, 0);// 0(存在), -1(不存在)
            if(isExistSubDir == -1){
                char mkdir[512];
                sprintf(mkdir, "mkdir %s", vsubsubDir);
                system(mkdir);
            }else{
                char mkdir[512], rm[512];
                sprintf(mkdir, "mkdir %s", vsubsubDir);
                sprintf(rm, "rm -r %s", vsubsubDir);

                system(rm);//先delete
                system(mkdir);// 再mkdir
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
    
    // write tecplot head file
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
    // write node and variable
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
                    sqrt(vert->Velocity.X*vert->Velocity.X+vert->Velocity.Y*vert->Velocity.Y), // Velocity
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
// 后处理文件，针对 TaylorGreen 问题，计算误差, 2021.11.08
/*------------------------------------------------------------------------------------------------*/
// 压强解析解
double PressureAnalyticalSolution_TaylorGreen(double xc, double yc)
{
    double rho = 1;
    return (10 + rho / 4.0 *( cos(2*PI*xc) + cos(2*PI*yc) ) );
}


// 速度解析解
void VelocityAnalyticalSolution_TaylorGreen(double x, double y, double *vx, double *vy)
{
    *vx =   sin(PI*x)*cos(PI*y);
    *vy = - cos(PI*x)*sin(PI*y);
}


// 针对TaylorGreen问题，计算速度，压强的误差
void TaylorGreenL2Norm()
{
    int         v = 0, c = 0, vSize = MeshObj.VertsArrLen, cSize = MeshObj.TrgsArrLen;
    double      x = 0, y = 0, xc = 0, yc = 0;
    double      VelL2 = 0, PreL2 = 0;
    double      vx_exact, vy_exact, vx, vy, p_exact, p;

    /////// 节点量
    for  (v = 0; v < vSize; v++)
    {
        x = MeshObj.Vertices[v].Pos.X;
        y = MeshObj.Vertices[v].Pos.Y;
        
        // 计算速度解析解
        VelocityAnalyticalSolution_TaylorGreen(x, y, &vx_exact, &vy_exact);
        
        vx = MeshObj.Vertices[v].Velocity.X;
        vy = MeshObj.Vertices[v].Velocity.Y;

        VelL2 += (vx - vx_exact)*(vx - vx_exact) + (vy - vy_exact)*(vy - vy_exact);
    }
    VelL2 /= sqrt(vSize);


    /////// 单元量
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


// 针对TaylorGreen问题，计算速度，压强的 Volume-Weighted L1误差
void VolumeWeightedL1Norm()
{
    int         v = 0, c = 0, vSize = MeshObj.VertsArrLen, cSize = MeshObj.TrgsArrLen;
    double      x = 0, y = 0, xc = 0, yc = 0;
    double      VelL1 = 0, PreL1 = 0, TotalNodalVolume = 0, TotalElementVolume = 0;
    double      vx_exact, vy_exact, vx, vy, p_exact, p, elem_vol, *node_vol = NULL;

    node_vol = (double *) calloc(vSize, sizeof(double));
    /////// 单元量, 计算节点控制体的体积
    for(c = 0; c < cSize; c++)
    {
        xc = (MeshObj.Trgs[c].Vertices[0]->Pos.X + MeshObj.Trgs[c].Vertices[1]->Pos.X + MeshObj.Trgs[c].Vertices[2]->Pos.X) / 3.0;
        yc = (MeshObj.Trgs[c].Vertices[0]->Pos.Y + MeshObj.Trgs[c].Vertices[1]->Pos.Y + MeshObj.Trgs[c].Vertices[2]->Pos.Y) / 3.0;
        
        // printf("trg[%d], {(%.2f, %0.2f), (%.2f, %0.2f), (%.2f, %0.2f)}\n", c, MeshObj.Trgs[c].Vertices[0]->Pos.X, MeshObj.Trgs[c].Vertices[0]->Pos.Y, 
        // MeshObj.Trgs[c].Vertices[1]->Pos.X, MeshObj.Trgs[c].Vertices[1]->Pos.Y, MeshObj.Trgs[c].Vertices[2]->Pos.X, MeshObj.Trgs[c].Vertices[2]->Pos.Y);

        p_exact = PressureAnalyticalSolution_TaylorGreen(xc, yc);
        p = MeshObj.Trgs[c].Pressure;
        elem_vol = MeshObj.Trgs[c].Area;
        
        // 计算节点控制体的体积
        node_vol[MeshObj.Trgs[c].Vertices[0]->Index] += elem_vol / 3.0;
        node_vol[MeshObj.Trgs[c].Vertices[1]->Index] += elem_vol / 3.0;
        node_vol[MeshObj.Trgs[c].Vertices[2]->Index] += elem_vol / 3.0;

        if (MeshObj.Trgs[c].MaterialId > 0){
            // 计算压强的 Volume-Weighted L1误差
            TotalElementVolume += elem_vol;
            PreL1 += fabs(p_exact - p)*elem_vol;
        }
    }
    PreL1 /= TotalElementVolume;

    /////// 节点量
    for  (v = 0; v < vSize; v++)
    {
        x = MeshObj.Vertices[v].Pos.X;
        y = MeshObj.Vertices[v].Pos.Y;
        
        // 计算速度解析解
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

    /////// 节点量
    sprintf(vfileName, "%s.txt", "Velocity_x_BottomEdge");
    fp = fopen(vfileName, "w");

    fprintf(fp, "x\tvelocity_x\texact\n");
    for  (v = 0; v < vSize; v++)
    {
        if (!MeshObj.Vertices[v].IsOnDownBoundary) continue;
        
        x = MeshObj.Vertices[v].Pos.X;
        y = MeshObj.Vertices[v].Pos.Y;
        
        // 计算速度解析解
        VelocityAnalyticalSolution_TaylorGreen(x, y, &vx_exact, &vy_exact);
        
        vx = MeshObj.Vertices[v].Velocity.X;
        vy = MeshObj.Vertices[v].Velocity.Y;

        fprintf(fp, "%e\t%e\t%e\n", x, vx, vx_exact);
    }
    ////
    fclose(fp);

       
    /////// 单元量
    // 生成顶点相邻的三角形编号数组
    int *vnbtrgmark     = (int *) calloc(vSize, sizeof(int));
    int (*vnbtrg)[25]   = (int(*)[25]) calloc(vSize*25, sizeof(int));
    // int (*vnbtrg)[25]   = (int *) calloc(vSize*25, sizeof(int));
    int k, t;

	/// 遍历三角形单元
	for (c = 0; c < cSize; c++)
	{
		for (v = 0; v < 3; v++ )
		{
			k = MeshObj.Trgs[c].Vertices[v]->Index; // v节点编号
			t = vnbtrgmark[k];
			vnbtrg[k][t] = c; //存放相邻三角形单元编号
			vnbtrgmark[k]++;  //计算相邻三角形单元个数		
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
        
        // 计算压强解析解
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
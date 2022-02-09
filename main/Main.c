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
	// 全局变量初始化
	SetInitialOrDefaultValueOfGlobals();

	// 命令行获取参数
	SetParameter(argc, argv);
		
	// 导入计算模型
	ImportModel(modelInfoFileName);

	// 打印基本信息
	PrintInfo(); 
	PrintMeshScale();

	// 创建输出文件
	CreateOutPutFile();

	// 演化
	printf("Evolving: \r\n");
	Start_time = clock();
	Evolve();
	End_time = clock();
	TotTime = GetWallTime(Start_time, End_time);

	//////////////////////////////////////////////// 打印信息
	PrintTime();

	if(TaylorGreen){	
		TaylorGreenOutPutBottomEdgePhysics();
		VolumeWeightedL1Norm();
		TaylorGreenL2Norm();
	}	
}


void PrintHelpInfo()
{
	printf("\n");
	printf("-f  : Model file.\n");
	printf("-mf : Is there a matter flow method (default %d), 1: yes, 0: no.\n", HAVE_MF);
	printf("-mfm: Types of matter flow methods (default %d).\n", MF_Method);
	printf("      0: 中点的动量, 即 Delta P=（v1+v2）M/3;\n");
	printf("      1：中点 + Delta P （优化问题）;\n");
	printf("      2: Delta P （优化问题）.\n");
	printf("-tg : Is there Taylor-Green vortex problem (default %d), 1: yes, 0: no.\n", TaylorGreen);
	printf("-k  : Times of viscosity coefficient (default %.4f).\n", ViscCoeffTimes);
	printf("-cs : C_safe, time step factor (default %.4f).\n", timeStepSafeFactor);
	printf("-bn : How many pictures to output (default %d).\n", TotBitmapOutputs);
	printf("\n");
}

// 从命令行读取参数
void SetParameter(int argc, char *argv[])
{
	// printf("%s\n", argv[0]);//读取可执行程序（包括路径）
	/*读取参数*/
	int i = 1;
	while(i < argc)
	{	
		// printf("%s\n", argv[i]);
		if(!strcmp(argv[i], "-f"))
		{
			modelInfoFileName = argv[i + 1];
			// printf("%s %s\n", argv[i], argv[i + 1]);
		}
		if(!strcmp(argv[i], "-mf"))
		{
			HAVE_MF = atoi(argv[i + 1]);
			// printf("%s %d\n", argv[i], HAVE_MF);
		}
		if(!strcmp(argv[i], "-mfm"))
		{
			MF_Method = atoi(argv[i + 1]);
			// printf("%s %d\n", argv[i], HAVE_MF);
		}		
		if(!strcmp(argv[i], "-bn"))
		{
			TotBitmapOutputs = atoi(argv[i + 1]);
			// printf("%s %d\n", argv[i], TotBitmapOutputs);
		}	
		if(!strcmp(argv[i], "-k"))
		{
			ViscCoeffTimes = atof(argv[i + 1]);
			// printf("%s %d\n", argv[i], ViscCoeffTimes);
		}	
		if(!strcmp(argv[i], "-cs"))
		{
			timeStepSafeFactor = atof(argv[i + 1]);
			// printf("%s %d\n", argv[i], timeStepSafeFactor);
		}	
		if(!strcmp(argv[i], "-tg"))
		{
			TaylorGreen = atoi(argv[i + 1]);
			// printf("%s %d\n", argv[i], TaylorGreen);
		}	

		if(!strcmp(argv[i], "-help") || !strcmp(argv[i], "--help") || !strcmp(argv[i], "-h") || !strcmp(argv[i], "--h"))
		{
			printf("Option %s: \n", argv[i]);
			PrintHelpInfo();
			exit(-1);
		}		
		i += 2;
	}
}



void SetInitialOrDefaultValueOfGlobals()
{
	//////////////////////////////////////////// 格点和三角形数组初始化
	CreateVerticesArray(&MeshObj.Vertices, &MeshObj.DeadVertIdsArr);
	CreateTrgsArray(&MeshObj.Trgs, &MeshObj.DeadTrgIdsArr);
	//////////////////////////////////////////// 物理常数
	/// 气体普适常量
	ROfGas = 8.3149 * 1e-5;//8.3149*J.mol^-1.k^-1
	//////////////////////////////////////////// 运行时间控制参数初始化
	/// <summary>
	/// 记录系统演化的时间，用于续算.
	/// </summary>
	TimeEvolvedRecord = 0.0;
	/// <summary>
	/// 系统演化的时间数。（它等于时间除以某个间隔）
	/// </summary>
	TimeNumbers = 0;
	/// <summary>
	/// 用来控制数据输出的时间间隔
	/// </summary>
	TimeInterval = 1.0;
	/// <summary>
	/// 用来控制程序运行的结束时间
	/// </summary>
	TimeEnd = 1.0;
}

void CreateOutPutFile()
{
	PhysQuantStatisticsFile = MyOpenFile("Physical_Quantity_Statistics.txt", "w");
	/// 其它物理量统计文件
	{
		FILE * fp = PhysQuantStatisticsFile;
		fprintf(fp, "下表数据为各种物理量的统计");
		fprintf(fp, "\n\n");
		fprintf(fp, "     时间(us)  \t");
		/// 依次罗列各种希望输出的物理量
		{
			fprintf(fp, "  质量(g)      \t");
			fprintf(fp, "  动量x(g)     \t");
			fprintf(fp, "  动量y(g)     \t");
			fprintf(fp, "  能量(100KJ)  \t");
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
}

void WriteOutputData()
{
	/// 写统计型的物理量
	{
		int i;
		///////////////////////////////////////////////////////// 依次统计和输出各种希望输出的物理量统计
		double totMass = 0.0, totMomentumX = 0.0, totMomentumY = 0.0;
		double totEnergy = 0.0;
		for (i = 0; i < MeshObj.TrgsArrLen; i++)
		{
			Triangle trg = &MeshObj.Trgs[i];
			if (trg->IsDead) continue;
			///////////////////////////////////////
			totMass += trg->Mass;
		}
		for (i = 0; i < MeshObj.VertsArrLen; i++)
		{
			Vertex vert = &MeshObj.Vertices[i];
			if (vert->IsDead) continue;
			//////////////////////////////
			totMomentumX += vert->Mass * vert->Velocity.X;
		}
		for (i = 0; i < MeshObj.VertsArrLen; i++)
		{
			Vertex vert = &MeshObj.Vertices[i];
			if (vert->IsDead) continue;
			//////////////////////////////
			totMomentumY += vert->Mass * vert->Velocity.Y;
		}
		for (i = 0; i < MeshObj.TrgsArrLen; i++)
		{
			Triangle trg = &MeshObj.Trgs[i];
			if (trg->IsDead) continue;
			//////////////////////////////
			if (trg->MaterialId > 0)
			{
				totEnergy += trg->InternalEnergy;
			}
		}
		for (i = 0; i < MeshObj.VertsArrLen; i++)
		{
			Vertex vert = &MeshObj.Vertices[i];
			if (vert->IsDead) continue;
			//////////////////////////////
			totEnergy += 0.5 * vert->Mass * CalcLengthSquar(vert->Velocity);
			totEnergy += vert->Mass * vert->Pos.Y * MeshObj.GravityFactor;//加上重力势能
		}
		/// 然后逐个输出
		fprintf(PhysQuantStatisticsFile, "%15e", TimeEvolvedRecord);
		fprintf(PhysQuantStatisticsFile, "\t%15e", totMass);
		fprintf(PhysQuantStatisticsFile, "\t%15e", totMomentumX);
		fprintf(PhysQuantStatisticsFile, "\t%15e", totMomentumY);
		fprintf(PhysQuantStatisticsFile, "\t%15e", totEnergy);
		/////////////////////////////////////////////////////////
		fprintf(PhysQuantStatisticsFile, "\n");
		fflush(PhysQuantStatisticsFile);
	}
}

void ImportModel(char * modelFile)
{
	////////////////////////////////////////// 从建模文件导入物理条件和网格
	printf("Loading compute model information from file \"%s\" \r\n\r\n", modelFile);
	MeshObjectFromLoadModel(modelFile);
	printf("Compute model information loaded. \r\n\r\n");
}


void Evolve()
{
	clock_t t1, t2, temp_time;
	double deltT, Time;

	SetAllDependentVariablesOfTrgs();

	//////////////////////////////////////////////////////////////////////////////////////// 输出初始状态网格文件和画图
	printf("%3d: Iter = %6d, t = %-.6f\n", TimeCountForBitmap, Iter, TimeEvolvedRecord);
	WriteOutputData();
	WriteTecPlotCellDataFromNormalGrid();

	////////////////////////////////////////// 数据准备 ///////////////////////////////////////////
	/// 流动步的强度量（力）
	VelocityPosEdgeCondition();
	/// 设置重力 
	SetGravity();
	/// 计算格点力 
	CalculateVertexForce();
	/// 设置边界受力 
	ForceEdgeCondition();
	/// 确定时间步长  
	deltT = DetermineDeltT();

	////////////////////////////////////////// T = T + \Delta T ///////////////////////////////////////////
	while (TimeEvolvedRecord <= TimeEnd)
	{
		Iter++;

	    /// 流动步演化  
		DynamicEvolve(deltT);
	    /// 设置因变量 
	    SetAllDependentVariablesOfTrgs();
	    /// 物质流动演化  
		if(HAVE_MF){
			t1 = clock();
			MatterFlowEvolve(deltT);
			t2 = clock();
			MatterFlowTime += GetWallTime(t1, t2);
		}

		//////////////////////////////////////////// 为下一时间步做准备
	    /// 流动步的强度量（力）
		VelocityPosEdgeCondition();
		/// 设置重力 
		SetGravity();
		/// 计算格点力 
	    CalculateVertexForce();
		/// 设置边界受力 
	    ForceEdgeCondition();
	    /// 确定时间步长  
		deltT = DetermineDeltT();

	    //////////////////////////////////////////// 后处理
	    TimeEvolvedRecord += deltT;
		printf("\r  Iter = %6d, t = %10.6e, deltT = %10.6e", Iter, TimeEvolvedRecord, deltT);
	    /////////////////////////////////////////// 输出文件
	    if (TimeEvolvedRecord / (TimeEnd / TotBitmapOutputs) > TimeCountForBitmap + 1)
	    {
	        TimeCountForBitmap++;
	        //////////////////////////////// 输出
			temp_time = clock();
			Time = GetWallTime(Start_time, temp_time);
			printf("\r%3d: Iter = %6d, t = %10.6e, deltT = %10.6e, TimeUsed=%.2fs\n", TimeCountForBitmap, Iter, TimeEvolvedRecord, deltT, Time);
			WriteTecPlotCellDataFromNormalGrid();
	    }
	    if (TimeEvolvedRecord / TimeInterval > TimeNumbers + 1)
	    {
	        TimeNumbers++;
	        //////////////////////////////////////////// 写物理量到数据文件
	        WriteOutputData();
	    }
	}
	
}








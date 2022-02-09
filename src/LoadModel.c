/*! \file  LoadModel.c
 *
 *  \brief Define some load model functions
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"


/// <summary>
/// 读取 Ansys 网格文件名，包括格点文件名和单元文件名
/// </summary>
/// <param name="sr"></param>
void LoadMeshFileNames(FILE * fr, char nodesFileName[], char cellsFileName[])
{
	char line[1001], words[3][1001];
	/////////////////////////////////////////////
    nodesFileName[0] = '\0';
    cellsFileName[0] = '\0';
	/////////////////////////////////////////////
    while (true)
    {
		fgets(line, 1000, fr);
		if(feof(fr))
        {
            ERR("读取LoadMeshFileNames出错：文件格式不符合要求\n");
        }
        /////////////////////////////////////////
		words[1][0] = '\0';        
		sscanf(line, "%s%s%s", words[0], words[1], words[2]);
		if(strcmp(words[1], "=")==0)
		{
			if (strcmp(words[0], "NodesFileName")==0)
			{
				if(words[2][0] == '\"')
					strcpy(nodesFileName, words[2]+1);
				else
					strcpy(nodesFileName, words[2]);
				///////////////////////////////
				{
					size_t len = strlen(nodesFileName);
					if(nodesFileName[len-1] == '\"')
					{
						nodesFileName[len-1] = '\0';
					}
				}
			}
			else if (strcmp(words[0], "CellsFileName")==0)
			{
				if(words[2][0] == '\"')
					strcpy(cellsFileName, words[2]+1);
				else
					strcpy(cellsFileName, words[2]);
				///////////////////////////////
				{
					size_t len = strlen(cellsFileName);
					if(cellsFileName[len-1] == '\"')
					{
						cellsFileName[len-1] = '\0';
					}
				}
				///////////////////////////////
				break;
			}
		}
    }
}
/// <summary>
/// 读取计算区域边界条件
/// </summary>
/// <param name="sr"></param>
void LoadBoxBoundaryCondition(FILE * fr, BoxBoundaryCondition * topBottomCondition, BoxBoundaryCondition * leftRightCondition)
{
	char line[1001], words[3][101];
	////////////////////////////////////////////////////////////////////
	(*topBottomCondition) = BoxBoundaryCondition_Wall;
	(*leftRightCondition) = BoxBoundaryCondition_Wall;
	////////////////////////////////////////////////////////////////////
	if(SearchNextInFile(fr, ">> 计算区域边界条件") == false)
	{
		ERR("没有找到“>> 计算区域边界条件”");
	}
	////////////////////////////////////////////////////////////////////
    while (true)
    {
        fgets(line, 1000, fr);
		if(feof(fr))
        {
            ERR("读取“计算区域边界条件”出错：文件格式不符合要求\n");
        }
        /////////////////////////////////////////
        words[1][0] = '\0';
		sscanf(line, "%s%s%s", words[0], words[1], words[2]);
		if(strcmp(words[1], "=")==0)
		{
			if (strcmp(words[0], "TopBottomBoundaryCondition") == 0)
			{
				if (strcmp(words[2], "固壁") == 0)
				{
					(*topBottomCondition) = BoxBoundaryCondition_Wall;
				}
				else if (strcmp(words[2], "周期") == 0)
				{
					(*topBottomCondition) = BoxBoundaryCondition_Cycle;
				}
				else
				{
					ERR("读取“计算区域边界条件”出错：文件格式不符合要求\n");
				}
			}
			else if (strcmp(words[0], "LeftRightBoundaryCondition") == 0)
			{
				if (strcmp(words[2], "固壁") == 0)
				{
					(*leftRightCondition) = BoxBoundaryCondition_Wall;
				}
				else if (strcmp(words[2], "周期") == 0)
				{
					(*leftRightCondition) = BoxBoundaryCondition_Cycle;
				}
				else
				{
					ERR("读取“计算区域边界条件”出错：文件格式不符合要求\n");
				}
				break;
			}
		}
    }
}
/// <summary>
/// 读取重力加速度
/// </summary>
/// <param name="?"></param>
/// <param name="meshObject"></param>
double ReadGravityAcceleration(FILE * fr)
{
	char line[1001], words[3][101];
	////////////////////////////////////////////////////////////////////
    double gravityAcc = 0.0;
	////////////////////////////////////////////////////////////////////
	if(SearchNextInFile(fr, ">> 重力加速度") == false)
	{
		ERR("没有找到“>> 重力加速度”");
	}
	////////////////////////////////////////////////////////////////////
    while (true)
    {
        fgets(line, 1000, fr);
		if(feof(fr))
        {
            ERR("读取“重力加速度”出错：文件格式不符合要求\n");
        }
        /////////////////////////////////////////
        words[1][0] = '\0';
		sscanf(line, "%s%s%s", words[0], words[1], words[2]);
		if(strcmp(words[1], "=")==0 && strcmp(words[0], "GravityAcceleration")==0)
		{
			gravityAcc = StringToDouble(words[2]);
			break;
		}
    }
    return (gravityAcc);
}

/// <summary>
/// 读取 数据输出的时间间隔 和 程序结束时间
/// </summary>
/// <param name="sr"></param>
void LoadTimeInterval(FILE * fr, double * timeInterval, double * timeEnd)
{
	char line[1001], words[3][101];
	////////////////////////////////////////////////////////////////////
    (*timeInterval) = 0.0;
    (*timeEnd) = 0.0;
	////////////////////////////////////////////////////////////////////
	if(SearchNextInFile(fr, ">> 数据输出的时间间隔") == false)
	{
		ERR("没有找到“>> 数据输出的时间间隔”");
	}
	////////////////////////////////////////////////////////////////////
    while (true)
    {
		if(feof(fr))
        {
            ERR("读取“数据输出的时间间隔”出错：文件格式不符合要求\n");
        }
        fgets(line, 1000, fr);
        /////////////////////////////////////////
        words[1][0] = '\0';
		sscanf(line, "%s%s%s", words[0], words[1], words[2]);
		if(strcmp(words[1], "=")==0)
		{
			if(strcmp(words[0], "TimeIntervalForOutput")==0)
			{
				(*timeInterval) = StringToDouble(words[2]);
			}
			else if (strcmp(words[0], "TimeEnd")==0)
			{
				(*timeEnd) = StringToDouble(words[2]);
				break;
			}
		}
    }
}

/// <summary>
/// 从输入文件 modelFile 读取计算模型。
/// </summary>
/// <param name="modelFile">建模文件</param>
/// <returns></returns>
void MeshObjectFromLoadModel(char * modelFile)
{
	FILE * fr;
    char nodesFileName[101];
    char cellsFileName[101];
    ////////////////////////////////////////////////////////// 1)打开文件
    fr = MyOpenFile(modelFile, "r");
    ////////////////////////////////////////////////////////// 2) 读取材料列表
    LoadMaterialListFromFile(fr);
    ////////////////////////////////////////////////////////// 3) 读取网格建模
    LoadMeshFileNames(fr, nodesFileName, cellsFileName);
    SetMeshFromMeshFiles(nodesFileName, cellsFileName);
    //////////////////////////////////////////////////////////// 4) 读取计算区域边界条件、重力加速度
    LoadBoxBoundaryCondition(fr, &MeshObj.TopBottomBoundaryCondition, &MeshObj.LeftRightBoundaryCondition);
    MeshObj.GravityFactor = ReadGravityAcceleration(fr);
    //////////////////////////////////////////////////////////// 5) 重设所有因变量
    SetAllDependentVariablesOfTrgs();
    ////////////////////////////////////////////////////////////    趁此机会根据初始时刻的最大物质声速设定空腔的模量（不采用空腔模量动态调整的做法了）20181029
	{
		double maxSoundVelocity = 0.0;
		int i;
		for (i = 0; i <MeshObj.TrgsArrLen; i++)
		{
			Triangle trg =&MeshObj.Trgs[i];
			if (trg->IsDead) continue;
			//////////////////////////////
			if (trg->SoundVelocitySum > maxSoundVelocity)
			{
				maxSoundVelocity = trg->SoundVelocitySum;
			}
		}
		MatParasList[0].Kai = MatParasList[0].NormalDensity * maxSoundVelocity * maxSoundVelocity;
	}
    //////////////////////////////////////////////////////////// 6) 读取数据输出的时间间隔
    LoadTimeInterval(fr, &TimeInterval, &TimeEnd);
    
    fclose(fr);
}
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


/////////////////////////////////////////////////////////////////////////////////////////// 全局物理常数
/// 气体普适常量
double ROfGas = 8.3149 * 1.0e-5; // 8.3149*J.mol^-1.k^-1
//////////////////////////////////////////////////////////////////
/// Mesh structure
struct meshObject MeshObj;
/// <summary>
/// 第0个元素——无意义，第1个元素——绝缘体类型的墙壁材料，第2个元素——理想导体类型的墙壁材料。
/// </summary>
struct materialParameters MatWall[3];
/// 物质种类集合
struct materialParameters MatParasList[MAX_MAT_TYPE];
int MatParasListLen;
//////////////////////////////////////////////////////////////////
/// <summary>
/// 记录系统演化的时间，用于续算.
/// </summary>
double TimeEvolvedRecord = 0.0;
/// <summary>
/// 系统演化的时间数。（它等于时间除以某个间隔）
/// </summary>
int TimeNumbers = 0;
/// <summary>
/// 图片数（它等于时间除以图片输出间隔）
/// </summary>
int TimeCountForBitmap = 0;
/// <summary>
/// 用来控制数据输出的时间间隔
/// </summary>
double TimeInterval = 1.0;
/// <summary>
/// 用来控制程序运行的结束时间
/// </summary>
double TimeEnd = 1.0;
/// <summary>
/// 输出的图片数目
/// </summary>
int TotBitmapOutputs = 100;
int Iter = 0;
//////////////////////////////////////////////////////////////////
FILE * PhysQuantStatisticsFile;
char * modelInfoFileName = "./data/0.05/Compute_Model_For_TriAngels.txt";


//////////////////////////////// 物质流动 //////////////////////////////////
int HAVE_MF     = 1;     // 物质流动开关, 1: ON, 0: OFF
int MF_Method   = 1;     // 0: 中点的动量, 即 Delta P=（v1+v2）M/3； 1：中点 + Delta P （优化问题）, 缺省值;  2: Delta P （优化问题）
int PistonFlag  = 0;     // 活塞问题标记
int TaylorGreen = 1;     // 是否为Taylor-Green vortex Problem, 能量方程需要增加源项



// 时间
double TotTime = 0; // 总时间
double MatterFlowTime = 0;
clock_t Start_time = 0;
clock_t End_time = 0;


bool CellDirFirstCreateFlag = true;
bool VertexDirFirstCreateFlag = true;


// 粘性系数的倍数（可调）
double ViscCoeffTimes = 2; 
double timeStepSafeFactor = 0.5; // 0.5 时间步长因子


vec2D ZeroVec2D = { 0, 0 };
vec2D UnitVecR = { 1, 0 };
tensor2D UnitTensor = {1, 0, 0, 1};
tensor2D ZeroTensor = {0, 0, 0, 0};





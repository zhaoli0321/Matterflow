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
    /// <summary>
    /// 格点数组中的位置索引。
    /// </summary>
    int Index;
    /// <summary>
    /// 标记格点是否已是“失效”格点。
    /// 用于导入Ansys网格时做标记，留待之后再删除。
    /// </summary>
    bool IsDead;
    /// <summary>
    /// 用于标记格点是否处于上下底边
    /// </summary>
    bool IsOnUpDownBoundary;
    bool IsOnUpBoundary;    // zhaoli
    bool IsOnDownBoundary;  // zhaoli
    /// <summary>
    /// 用于标记格点是否处于左右边界（网格镜像翻倍后，轴心处的格点也被定义为左右边界格点）
    /// </summary>
    bool IsOnLeftRightBoundary;
    bool IsOnLeftBoundary;  // zhaoli
    bool IsOnRightBoundary; // zhaoli
    /// <summary>
    /// 格点质量
    /// </summary>
    double Mass;
    /// <summary>
    /// 格点速度
    /// </summary>
    vec2D Velocity;
    /// <summary>
    /// 格点位置
    /// </summary>
    vec2D Pos;
    ////////////////////////////////// 临时辅助变量
    /// <summary>
    /// 格点受力
    /// </summary>
    vec2D Force;
    /// <summary>
    /// 格点受力中用于物质流动计算的部分
    /// </summary>
    vec2D ForceForMatterFlow;
    /// <summary>
    /// 一个时间步的位移。
    /// 用于辅助计算内能增量。
    /// </summary>
	vec2D DeltPos;
} *Vertex;


typedef struct triangle
{
    /// <summary>
    /// 三角形标号，用于调试程序。
    /// </summary>
    int Index;
    /// <summary>
    /// 标记三角形是否已是“失效”三角形。
    /// 用于Merge网格重分时暂时将该删除的三角形做标记，留待之后再删除。
    /// </summary>
    bool IsDead;
    /// <summary>
    /// 三角形的三个顶点（约定按逆时针方向排布）
    /// </summary>
    Vertex Vertices[3];
    /// <summary>
    /// 三个邻域三角形
    /// </summary>
    struct triangle * NeighbourTrgs[3];
    /// <summary>
    /// 与每条边对应的物质流出加速度
    /// </summary>
    double FlowAcc[3];
    /// <summary>
    /// 与每条边对应的补偿速度（物质流出速度）
    /// </summary>
    double FlowVelocity[3];
    /// <summary>
    /// 三角形的质量
    /// </summary>
    double Mass;
    /// <summary>
    /// 材料号（在材料列表中的位置——不是用户输入的材料编号）
    /// </summary>
    int MaterialId;
    /// <summary>
    /// 密度（中间辅助变量）
    /// </summary>
    double Density;
    /// <summary>
    /// 压强（中间辅助变量）
    /// </summary>
	double Pressure;
    /// <summary>
    /// 单元中的总模量确定的声速。（用于确定时间步长、物质流动法中确定衰减系数）
    /// </summary>
    double SoundVelocitySum;
    /// <summary>
    /// 动态人工粘性系数（中间辅助变量）
    /// </summary>
    double ViscCoeff;
    /// <summary>
    /// 粘性应力的 RZ 分量（中间辅助变量）
    /// </summary>
    /// 标量粘性应力  
    double Viscos;  
    /// <summary>
    /// 流体内应力之和。包括：压强、弹性应力、粘性应力。但不包含磁场应力张量
    /// （中间辅助变量）
    /// </summary>
	tensor2D TotalFluidStress;
    /// <summary>
    /// 三角形面积（中间辅助变量）
    /// </summary>
    double Area;
    /// <summary>
    /// 该三角单元三条边对应的三条高（中间辅助变量）
    /// </summary>
    double Heights[3];
    /// <summary>
    /// 该三角单元的最小高（中间辅助变量）
    /// </summary>
    double MinHeight;
    /// <summary>
    /// 格点坐标（中间辅助变量）
    /// </summary>
	vec2D CycledPoses[3];
    ////////////////////////////////////////////////////// 物质单元专用变量
    /// <summary>
    /// 单元中的流体内能
    /// </summary>
    double InternalEnergy;
    /// <summary>
    /// 温度（中间辅助变量）
    /// </summary>
    double Temperature;
} *Triangle;


typedef enum boxBoundaryCondition
{
	BoxBoundaryCondition_Wall,
	BoxBoundaryCondition_Cycle
} BoxBoundaryCondition;

typedef struct meshObject
{
    /////////////////////////////////////////////////// 网格
    /// <summary>
    /// 格点集合
    /// </summary>
	struct vertex * Vertices;
    int VertsArrLen;
    /// <summary>
    /// 失效格点集合
    /// </summary>
	int * DeadVertIdsArr;
    int DeadVertIdsArrLen;
    /// <summary>
    /// 三角形集合
    /// </summary>
	struct triangle * Trgs;
    int TrgsArrLen;
    /// <summary>
    /// 失效三角形集合
    /// </summary>
	int * DeadTrgIdsArr;
    int DeadTrgIdsArrLen;
    /////////////////////////////////////////////////// 计算区域边界格点条件（对于zr程序，左右必须是固壁边界条件）
    BoxBoundaryCondition TopBottomBoundaryCondition;
	BoxBoundaryCondition LeftRightBoundaryCondition;
    /////////////////////////////////////////////////// 全局条件
    /// <summary>
    /// 重力加速度
    /// </summary>
    double GravityFactor;
    //////////////////////////////// 计算区域大小
    double Width;
    double Height;
} *MeshObject;

typedef enum eOSTypeEnum { 
	EOSTypeEnum_Linear, 
	EOSTypeEnum_IdealGas
} EOSTypeEnum;

typedef struct materialParameters
{
    //////////////////////////////////////////////////////////////////////////////////////////// 材料属性
    int MaterialNumber;//对应用户输入文件中的材料号。
    //////////////////////////////////////////////////////////// 材料名称
    char MaterialName[101];
    //////////////////////////////////////////////////////////// EOS 参数
    EOSTypeEnum EOSType;     //材料类型
    /// 线弹性流体参数
    double NormalDensity;    //室温室压密度(共享TableEOS中的NormalDensity变量)
    double Kai;              //体模量
    double SpecificHeat;     //比热
    /// 多方气体参数
    double Gamma;            //多方指数
    double MolMass;          //摩尔质量
    //////////////////////////////////////////////////////////// 粘性参数
    double MotionViscos;
} *Material;

typedef struct nodeFromFile
{
	/// 格点坐标
	vec2D Pos;
	/// 格点速度
	vec2D Vel;
	/// 格点的邻接三角形
	int TrgIDs[MAX_NB_TRGS];
	int TrgIDsLen;
} *NodeFromFile;

typedef struct triangleFromFile
{
	/// 材料号
	int MaterialId;
    /// 压强
    double Pressure;
    /// 温度
    double Temperature;
	/// 三角形的邻接格点
	int NodeIDs[3];
	/// 三角形的邻接三角形
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

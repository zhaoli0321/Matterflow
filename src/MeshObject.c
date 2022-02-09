/*! \file  MeshObject.c
 *
 *  \brief Define some kernel functions
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"

void ERR(char * prnt)
{
    printf(prnt);
    getchar();
    exit(0);
}

double GetWallTime(clock_t start_time, clock_t end_time)
{
	return ( (double) (end_time - start_time) / CLOCKS_PER_SEC );
}

/*
	min{ |p| }
	a<p, p> + <b, p> + c = 0
 */	
vec2D CalculateMoveMomentumMatterFlow2(double a, vec2D b, double c)
{
	vec2D x_O;// x_O为圆心坐标
	double r, d; // r为圆的半径，d为点lambda到圆心的距离
	
	if(a < 0){
		printf("error a = %e < 0\n", a);
	}
	x_O = CValueMultVec2D(- 1.0 / (2*a), b);
	d = CalcLength(x_O);
	/////////////Radius
	double rr = CalcLengthSquar(b) / (4*a*a) - c / a;
	if(rr < 0){
		// printf("\nerror radius = %e < 0\n", rr);		
		return (x_O);
	}else if(rr == 0)
	{
		return (x_O);
	}else{
		r = sqrt(rr);
	} 

	if(d > 0)
	{
		return ( CValueMultVec2D(1 - r / d, x_O) );
	}else{
		// printf("error distance = %e <= 0\n", d);
		vec2D rn = {r, 0};
		// vec2D rn = {0.5 * sqrt(2) * r, 0.5 * sqrt(2) * r};
		return (rn);
	}
}

/*
	min{ |lambda - p| }
	a<p, p> + <b, p> + c = 0
 */
vec2D CalculateMoveMomentumMatterFlow1(double a, vec2D b, double c, vec2D lambda)
{
	vec2D x_O;// x_O为圆心坐标
	double r, d, ratio; // r为圆的半径，d为点lambda到圆心的距离
	
	if(a < 0){
		printf("error a = %e < 0\n", a);
	}
	x_O = CValueMultVec2D(- 1.0 / (2*a), b);
	d = CalcLength( Vec2DSub(lambda, x_O) );
	double rr = CalcLengthSquar(b) / (4*a*a) - c / a;
	if(rr <= 0)
	{
		// printf("\nerror radius = %e < 0\n", rr);
		// printf("a = %e, b^2 = %e, c = %e\n", a, CalcLengthSquar(b), c);
		// r = d;
		return (x_O);
	}else{
		r = sqrt(rr);
	} 

	if(d > 0)
	{
		ratio = r / d;
	}else{
		// printf("error distance = %e <= 0\n", d);
		vec2D rn = {r, 0};
		return ( Vec2DAdd(x_O, rn) );
	}
	/////////////
	return ( Vec2DAdd(CValueMultVec2D(ratio, lambda), CValueMultVec2D(1-ratio, x_O)) );
}


////// 二元二次函数
double func(vec2D x, double a, vec2D b, double c)
{
	return (a * Dot(x, x) + Dot(b, x) + c);
}



/// <summary>
/// 在Vertices的整个满数组被创建时，就把index初始化好，以后不再需要改变
/// </summary>
/// <param name="arrLen"></param>
/// <returns></returns>
void CreateVerticesArray(struct vertex * * vertsArrOut, int * * deadVertIdsArrOut)
{
	struct vertex * vertsArr = (struct vertex *)calloc(MAX_VERTS, sizeof(struct vertex));
	int * deadVertIdsArr = (int *)malloc(sizeof(int)*MAX_VERTS);
	int i;
	for (i = 0; i < MAX_VERTS; i++)
	{
		vertsArr[i].Index = i; 
		//MakeLockTagString(i, 10, vertsArr[i].LockTag);
	}
	(*vertsArrOut) = vertsArr;
	(*deadVertIdsArrOut) = deadVertIdsArr;
}

//////////////////////////////////////////////////////////////////////////////////////
/// <summary>
/// 在Trgs的整个满数组被创建时，就把index初始化好，以后不再需要改变
/// </summary>
/// <param name="arrLen"></param>
/// <returns></returns>
void CreateTrgsArray(struct triangle * * trgsArrOut, int * * deadTrgIdsArrOut)
{
	struct triangle * trgsArr = (struct triangle *)calloc(MAX_VERTS * 2, sizeof(struct triangle));
	int * deadTrgIdsArr = (int *)malloc(sizeof(int)*MAX_VERTS * 2);
	int i;
	for (i = 0; i < MAX_VERTS * 2; i++)
	{
		trgsArr[i].Index = i;
	}
	(*trgsArrOut) = trgsArr;
	(*deadTrgIdsArrOut) = deadTrgIdsArr;
}

/// <summary>
/// 
/// </summary>
/// <param name="vert"></param>
void SetDeadVertex(Vertex vert)
{
	vert->IsDead = true; // 设置失效标记
	MeshObj.DeadVertIdsArr[MeshObj.DeadVertIdsArrLen++] = vert->Index;
}

/// <summary>
/// 
/// </summary>
/// <param name="trg"></param>
void SetDeadTriangle(Triangle trg)
{
	trg->IsDead = true; // 设置失效标记
	MeshObj.DeadTrgIdsArr[MeshObj.DeadTrgIdsArrLen++] = trg->Index;
}

/// <summary>
/// 对一小群格点坐标（认为这些格点是相互邻近的）进行跨边界周期处理
/// </summary>
/// <param name="posesInit"></param>
/// <returns></returns>
// void PosesAfterCycleDealing(vec2D posesInit[], int totVertices, vec2D posesOut[])
// {
// 	int i;
// 	/////////////////////////////////////////////////////////////// 如果格点之间的坐标之差较大，则说明有的格点来自与跨界三角形，因此需要对来
// 	/////////////////////////////////////////////////////////////// 自于跨界三角形的格点的坐标进行改变
// 	for (i = 0; i < totVertices; i++)
// 	{
// 		posesOut[i] = posesInit[i];
// 	}
// 	///////////////////////////////////// 先对r坐标进行处理。
// 	/// 先前为了为了适应z-r坐标系，不允许出现r小于0的情形，
// 	/// 所以约定r坐标的周期处理中将坐标加上Width（以前是减去Width）。
// 	/// 现在，程序中采用了镜像翻倍的坐标，就允许r小于0，就还是可以约定减去width*2
// 	{
// 		int maxID;//最大r坐标顶点的编号
// 		double maxR;//最大r坐标
// 		maxR = posesOut[0].X;
// 		maxID = 0;
// 		for (i = 1; i < totVertices; i++)
// 		{
// 			if (posesOut[i].X > maxR)
// 			{
// 				maxR = posesOut[i].X;
// 				maxID = i;
// 			}
// 		}
// 		///然后将与最大的x坐标相差太大的顶点坐标加上Width
// 		for (i = 1; i < totVertices; i++)
// 		{
// 			if (maxR - posesOut[(maxID + i) % totVertices].X > MeshObj.Width * WidthFactor) posesOut[(maxID + i) % totVertices].X += MeshObj.Width;
// 		}
// 	}
// 	///////////////////////////////////// 然后对y坐标进行同样的处理
// 	{
// 		int minYID;//最小y坐标顶点的编号
// 		double minY;//最小y坐标
// 		minY = posesOut[0].Y;
// 		minYID = 0;
// 		for (i = 1; i < totVertices; i++)
// 		{
// 			if (posesOut[i].Y < minY)
// 			{
// 				minY = posesOut[i].Y;
// 				minYID = i;
// 			}
// 		}
// 		///然后将与最小的y坐标相差太大的顶点坐标减去Height
// 		for (i = 1; i < totVertices; i++)
// 		{
// 			if (posesOut[(minYID + i) % totVertices].Y - minY > MeshObj.Height * HeightFactor) posesOut[(minYID + i) % totVertices].Y -= MeshObj.Height;
// 		}
// 	}
// }

/// <summary>
/// 用周期条件将格点的坐标限制在(-@idth,Width)(0,Height)之间。
/// </summary>
/// <param name="pos"></param>
/// <returns></returns>
// vec2D RestrainPosWithCycleCondition(vec2D pos)
// {
// 	if (pos.X > MeshObj.Width)
// 	{
// 		pos.X -= MeshObj.Width;
// 	}
// 	else if (pos.X < 0.0)
// 	{
// 		pos.X += MeshObj.Width;
// 	}
// 	///////////////////////////
// 	if (pos.Y > MeshObj.Height)
// 	{
// 		pos.Y -= MeshObj.Height;
// 	}
// 	else if (pos.Y < 0.0)
// 	{
// 		pos.Y += MeshObj.Height;
// 	}
// 	///////////////////////////
// 	return (pos);
// }


/// <summary>
/// 给每个Vertex设置重力
/// </summary>
void SetGravity()
{
	int VertsArrLen = MeshObj.VertsArrLen;
	struct vertex * Vertices = MeshObj.Vertices;
	int k;
	Vertex vert;
	///////////////////////////////////////
	for (k = 0; k < VertsArrLen; k++)
	{
		vert = &Vertices[k];
		if (vert->IsDead) continue;
		//////////////////////////////////////
		vert->Force = MakeVec2D(0, -MeshObj.GravityFactor * vert->Mass);
		vert->ForceForMatterFlow = ZeroVec2D;
	}
}


/// <summary>
/// 计算格点力
/// </summary>
void CalculateVertexForce()
{
	int VertsArrLen = MeshObj.VertsArrLen, TrgsArrLen = MeshObj.TrgsArrLen;
	Vertex Vertices = MeshObj.Vertices;
	Triangle Trgs = MeshObj.Trgs;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	int k, i;
	Triangle trg;
	Vertex vert;
	Vertex* trgVertices;
	tensor2D stress;
	vec2D side, f;
	///////////////////////////////////////
	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		if (trg->IsDead) continue;
		//////////////////////////////
		if (trg->MaterialId < 0) continue;
		///////////////////////////////////////////////////////////////////////////////////////////////// 如果是物质单元，要包含磁压
		{
			trgVertices = trg->Vertices;
			if (trg->MaterialId > 0)
			{
				////////////////////////////////////////////////////////////////// 计算用于物质流动的受力
				{						
					stress = CValueToTensor2D(trg->Pressure + trg->Viscos);					
					///////////////////////////////////// 把应力分到三个格点上
					for (i = 0; i < 3; i++)
					{
						/////////////////////////////// 三角形第i条边对应的圆台面
						side = Rotate270(Vec2DSub(trg->CycledPoses[(i + 1) % 3], trg->CycledPoses[i]));
						/////////////////////////////// 应力分到格点上
						/// 全部力
						f = Tensor2DMultVec(stress, side);
						Vec2DSelfAdd(&trgVertices[i]->ForceForMatterFlow, Vec2DMultCValue(f, 0.5));
						Vec2DSelfAdd(&trgVertices[(i + 1) % 3]->ForceForMatterFlow, Vec2DMultCValue(f, 0.5));
					}
				}
				////////////////////////////////////////////////////////////////// 计算完整受力
				{
					stress = trg->TotalFluidStress;
					///////////////////////////////////// 把应力分到三个格点上
					for (i = 0; i < 3; i++)
					{
						/////////////////////////////// 三角形第i条边对应的圆台面
						side = Rotate270(Vec2DSub(trg->CycledPoses[(i + 1) % 3], trg->CycledPoses[i]));
						/////////////////////////////// 应力分到格点上
						/// 全部力
						f = Tensor2DMultVec(stress, side);
						Vec2DSelfAdd(&trgVertices[i]->Force, Vec2DMultCValue(f, 0.5));
						Vec2DSelfAdd(&trgVertices[(i + 1) % 3]->Force, Vec2DMultCValue(f, 0.5));
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////// 如果是真空单元
			else
			{
				stress = trg->TotalFluidStress;
				///////////////////////////////////// 把应力分到三个格点上
				for (i = 0; i < 3; i++)
				{
					/////////////////////////////// 三角形第i条边对应的圆台面
					side = Rotate270(Vec2DSub(trg->CycledPoses[(i + 1) % 3], trg->CycledPoses[i]));
					/////////////////////////////// 应力分到格点上
					/// 全部力
					f = Tensor2DMultVec(stress, side);
					Vec2DSelfAdd(&trgVertices[i]->Force, Vec2DMultCValue(f, 0.5));
					Vec2DSelfAdd(&trgVertices[(i + 1) % 3]->Force, Vec2DMultCValue(f, 0.5));
				}
			}
		}
	}
}


/// <summary>
/// 位置速度边界条件
/// </summary>
/// <param name="time"></param>
void VelocityPosEdgeCondition()
{
	int VertsArrLen = MeshObj.VertsArrLen, TrgsArrLen = MeshObj.TrgsArrLen;
	Vertex Vertices = MeshObj.Vertices;
	Triangle Trgs = MeshObj.Trgs;
	BoxBoundaryCondition TopBottomBoundaryCondition = MeshObj.TopBottomBoundaryCondition;
	double Width  = MeshObj.Width;
	double Height = MeshObj.Height;
	/////////////////////////////////////////////////////////////////////////
	////////////// 边界格点满足的边界条件
	int k;
	Vertex vert;
	///////////////////////////////////////
	for (k = 0; k < VertsArrLen; k++)
	{
		vert = &Vertices[k];
		if (vert->IsDead) continue;
		/////////// 如果底边是固壁边界条件，则底边格点的纵向坐标为零
		if (vert->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall)
		{
			vert->Velocity.Y = 0;
			if (vert->IsOnUpBoundary)	vert->Pos.Y = Height;	
			if (vert->IsOnDownBoundary)	vert->Pos.Y = 0.0;							
		}
		/////////// 左右边界格点的横向坐标为零
		if (vert->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall)
		{
			vert->Velocity.X = 0;
			if (vert->IsOnRightBoundary)	vert->Pos.X = Width;				
			if (vert->IsOnLeftBoundary)		vert->Pos.X = 0.0;							
		}
	}

	////////////// 墙壁单元上的格点满足的速度条件
	int t;
	Triangle trg;
	///////////////////////////////////////
	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		if (trg->IsDead) continue;
		//////////////////////////////
		if (trg->MaterialId < 0)
		{
			for (t = 0; t < 3; t++)
			{
				vert = trg->Vertices[t];
				vert->Velocity = ZeroVec2D;
			}
		}
	}
}

/// <summary>
/// 受力边界条件
/// </summary>
/// <param name="time"></param>
void ForceEdgeCondition()
{
	int VertsArrLen = MeshObj.VertsArrLen, TrgsArrLen = MeshObj.TrgsArrLen;
	Vertex Vertices = MeshObj.Vertices;
	Triangle Trgs = MeshObj.Trgs;
	BoxBoundaryCondition TopBottomBoundaryCondition = MeshObj.TopBottomBoundaryCondition;
	
	////////////// 边界格点满足的边界条件
	int k;
	Vertex vert;
	///////////////////////////////////////

	for (k = 0; k < VertsArrLen; k++)
	{
		vert = &Vertices[k];
		if (vert->IsDead) continue;
		///////////// 底边格点的纵向受力为零
		if (vert->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall)
		{
			vert->Force.Y = 0;
			vert->ForceForMatterFlow.Y = 0;
		}
		/////////// 左右边界格点的横向受力为零
		if (vert->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall)
		{
			vert->Force.X = 0;
			vert->ForceForMatterFlow.X = 0;
		}
	}

	////////////// 墙壁单元上的格点满足的受力条件
	int t;
	Triangle trg;
	///////////////////////////////////////
	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		if (trg->IsDead) continue;
		//////////////////////////////
		if (trg->MaterialId < 0)
		{
			for (t = 0; t < 3; t++)
			{
				vert = trg->Vertices[t];
				vert->Force = ZeroVec2D;
			}
		}
	}
}

/// <summary>
/// 系统演化一个时间步
/// </summary>
/// <param name="deltT"></param>
void DynamicEvolve(double deltT)
{
	int VertsArrLen = MeshObj.VertsArrLen, TrgsArrLen = MeshObj.TrgsArrLen;
	Vertex Vertices = MeshObj.Vertices;
	Triangle Trgs = MeshObj.Trgs;
	/////////////////////////////////////////////////////// 更新坐标和速度
	int k;
	Vertex vert;
	vec2D deltVelocity;
	///////////////////////////////////////
	for (k = 0; k < VertsArrLen; k++)
	{
		vert = &Vertices[k];
		if (vert->IsDead) continue;
		//////////////////////////////////////
		{
			deltVelocity = Vec2DMultCValue(vert->Force, 1 / vert->Mass * deltT);
			vert->DeltPos = Vec2DMultCValue(Vec2DAdd(vert->Velocity, CValueMultVec2D(0.5, deltVelocity)), deltT);//趁此机会记录位移增量，用于下面计算内能增量
			Vec2DSelfAdd(&vert->Pos, vert->DeltPos);
			Vec2DSelfAdd(&vert->Velocity, deltVelocity);
		}			
	}

	/////////////////////////////////////////////////////// 更新流体内能（磁场应力对流体内能无贡献）
	int i;
	Triangle trg;
	Vertex *verts;
	vec2D side, f;
	vec2D * poses;
	tensor2D stress;
	vec2D fLats[3];
	double xc = 0, yc = 0, Se = 0;
	///////////////////////////////////////
	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		if (trg->IsDead) continue;
		//////////////////////////////
		if (trg->MaterialId > 0)
		{
			/////////////////////////////////////////////////////////////////////////////// 将本三角形的内应力转换为三个格点的受力
			poses = trg->CycledPoses;
			stress = trg->TotalFluidStress;
			fLats[0].X = 0;
			fLats[0].Y = 0;
			fLats[1].X = 0;
			fLats[1].Y = 0;
			fLats[2].X = 0;
			fLats[2].Y = 0;

			for (i = 0; i < 3; i++)
			{
				/////////////////////////////// 三角形第i条边对应的圆台面
				side = Rotate270(Vec2DSub(trg->CycledPoses[(i + 1) % 3], trg->CycledPoses[i]));
				/////////////////////////////// 应力分到格点上
				/// 全部力
				f = Tensor2DMultVec(stress, side);
				Vec2DSelfAdd(&fLats[i], Vec2DMultCValue(f, 0.5));
				Vec2DSelfAdd(&fLats[(i + 1) % 3], Vec2DMultCValue(f, 0.5));
			}
			/////////////////////////////////////////////////////////////////////////////// 计算内能增量(应力做功)
			{
				verts = trg->Vertices;
				trg->InternalEnergy += (-Dot(fLats[0], verts[0]->DeltPos) - Dot(fLats[1], verts[1]->DeltPos) - Dot(fLats[2], verts[2]->DeltPos));
			}
			/////////////////////////////////////////////////////////////////////////////// 计算内能增量（源项）2021.05.15
			if (TaylorGreen){			
				// Taylor Green problem, add the energy source term S(x,y):
				// Se(x,y) = 3*pi/8 [cos(pi*x)cos(3*pi*y) - cos(3*pi*x)cos(pi*y)]
				// where (x,y) is the center of gravity of the element.
				xc = (verts[0]->Pos.X + verts[1]->Pos.X + verts[2]->Pos.X) / 3.0;
				yc = (verts[0]->Pos.Y + verts[1]->Pos.Y + verts[2]->Pos.Y) / 3.0;
				
				// Se = 3*PI/8.0 *( cos(PI*xc)*cos(3*PI*yc) - cos(3*PI*xc)*cos(PI*yc) );
				Se = 3*PI/4.0 *( cos(PI*xc)*sin(PI*yc)*sin(2*PI*yc) - sin(PI*xc)*sin(2*PI*xc)*cos(PI*yc) );
				// printf("Se = %f(Se*Mass = %f), InternalEnergy = %f, Mass = %f\n", Se, Se*trg->Mass, trg->InternalEnergy, trg->Mass);
				trg->InternalEnergy += Se * deltT * trg->Area;
				// printf("Se = %f, Se * deltT * trg->Area = %e\n", Se, Se * deltT * trg->Area);//trg->Area
			}
		}
	}	
}

/// <summary>
/// “物质流动补偿”。
/// 方案二（20170505）：“假想中间格点造成的边弯曲”方案。
/// 仅由压强差来确定物质流动量。
/// 计入物质流动对比体积变化的贡献，进而对粘性大小的影响。
/// //////////////////////////////////////////////////////
/// 本函数计算其中的物质流动加速度
/// </summary>
/// <param name="halfDeltT"></param>
void CalculateMatterFlowAcc()
{
	int TrgsArrLen = MeshObj.TrgsArrLen;
	Triangle Trgs = MeshObj.Trgs;
	BoxBoundaryCondition TopBottomBoundaryCondition = MeshObj.TopBottomBoundaryCondition;
	/////////////////////////////////////////////////////////////////////////////////////////////
	int i, j, k, t;
	Triangle trg, trgNb;
	Vertex vert, vert1, vert2;
	int idNb;
	vec2D side, sn, a1, a2, force, forceNb;
	double a12, ac;

	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		if (trg->IsDead) continue;
		//////////////////////////////
		if (trg->MaterialId <= 0) continue;//真空和墙壁单元不计算物质流动
		//////////////////////////////////////////// 遍历三角形的每条边
		for (i = 0; i < 3; i++)
		{
			///////////////////////////////////////////////////////////////////////////////////// 从边格点属性判断边的属性 
			vert1 = trg->Vertices[i];
			vert2 = trg->Vertices[(i + 1) % 3];
			/// 上下底边之间在固壁边界条件下不发生物质流动补偿
			if (vert1->IsOnUpDownBoundary && vert2->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowAcc[i] = 0.0;
				continue;
			}
			/// 左右边界之间在固壁边界条件下不发生物质流动补偿
			if (vert1->IsOnLeftRightBoundary && vert2->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowAcc[i] = 0.0;
				continue;
			}
			///////////////////////////////////////////////////////////////////////////////////// 找到该条边的相关邻域三角形、顶点和边信息
			trgNb = trg->NeighbourTrgs[i];			
			if(trgNb == NullTriangle) continue;
			
			/// 找到邻域三角形相关的顶点号
			idNb = -1;
			for (j = 0; j < 3; j++)
			{
				if (trgNb->NeighbourTrgs[j] == trg)
				{
					idNb = j;
					break;
				}
			}
			///////////////////
			if (trgNb->MaterialId != trg->MaterialId)// 不同种材料之间不发生物质流动补偿
			{
				trg->FlowAcc[i] = 0.0;
				continue;
			}
			{
				///////////////////////////////////////////////////////////////////////////////////// 计算出边相关的信息备用
				side = Rotate90(Vec2DSub(trg->CycledPoses[i], trg->CycledPoses[(i + 1) % 3]));
				sn = UnitVec(side);
				///////////////////////////////////////////////////////////////////////////////////// 计算三角形边上的流出速率加速度
				//////////////////// 两个顶点的加速度和它们的平均加速度
				a1 = Vec2DDivideCValue(vert1->ForceForMatterFlow, vert1->Mass);
				a2 = Vec2DDivideCValue(vert2->ForceForMatterFlow, vert2->Mass);
				a12 = Dot(CValueMultVec2D(0.5, Vec2DAdd(a1, a2)), sn);//要去掉重力加速度
				//////////////////// 设想的中间格点的加速度						
				force 	= CValueMultVec2D((trg->Pressure + trg->Viscos) * 0.5, side);
				forceNb = CValueMultVec2D((trgNb->Pressure + trgNb->Viscos) * (-0.5), side);					
				//double ac = Dot(force + forceNb, sn) / ((trg->Mass + trgNb->Mass) / 3.0);
				ac = Dot(Vec2DAdd(force, forceNb), sn) / ((trg->Mass + trgNb->Mass) / 4.0);//谢龙提议中心格点的质量用“1/4”的三角形质量。20190921
				//////////////////// 
				trg->FlowAcc[i] = ac - a12;
			}
		}
	}
}


/// <summary>
/// “物质流动”。
/// </summary>
/// <param name="halfDeltT"></param>
void MatterFlowEvolve(double deltT)
{	
	int TrgsArrLen = MeshObj.TrgsArrLen;
	Triangle Trgs = MeshObj.Trgs;
	BoxBoundaryCondition TopBottomBoundaryCondition = MeshObj.TopBottomBoundaryCondition;
	/////////////////////////////////////////////////////////////////////////////////////////////
	double halfDeltT = deltT * 0.5; //由于每条边要被遍历两次，因此将时间减半

	CalculateMatterFlowAcc();

	/////////////////////////////////////////////////////////////////////////////// 一边循环，一边就修改了单元和格点的质量、能量、速度
	int i, j, k;
	Triangle trg = NULL, trgNb = NULL;
	Vertex vert1, vert2, vertOpp, vertNbOpp;
	double massFlowOut;
	int idNb = -1;
	double sideLength, moveDistance, factorDecrease, flowDensity;
	double internalEnergy, internalEnergyNb, energyFlowOut, initialKineticEnergy, deltMassVert, deltKineticEnergy;
	vec2D initialMomentumOpp, initialMomentumNbOpp, momentumFlowOut, newMomentumOpp, newMomentumNbOpp, deltMomentumVert;
	double a, c;
	vec2D b, lambda;

	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		if (trg->IsDead) continue;
		//////////////////////////////
		if (trg->MaterialId <= 0) continue;//真空和墙壁单元不计算物质流动
		//////////////////////////////////////////// 遍历三角形的每条边
		for (i = 0; i < 3; i++)
		{
			///////////////////////////////////////////////////////////////////////////////////// 找到该条边的相关邻域三角形、顶点和边信息
			trgNb = trg->NeighbourTrgs[i];			
			if(trgNb == NullTriangle) continue;

			/// 找到邻域三角形相关的顶点号
			for (idNb = -1, j = 0; j < 3; j++)
			{
				if (trgNb->NeighbourTrgs[j] == trg)
				{
					idNb = j;
					break;
				}
			}
			///////////////////
			if (trgNb->MaterialId != trg->MaterialId)// 不同种材料之间不发生物质流动补偿
			{
				trg->FlowVelocity[i] = 0;
				trgNb->FlowVelocity[idNb] = 0;
				continue;
			}
			/////////////////// 
			vert1 = trg->Vertices[i];
			vert2 = trg->Vertices[(i + 1) % 3];
			/// 上下底边之间在固壁边界条件下不发生物质流动补偿
			if (vert1->IsOnUpDownBoundary && vert2->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowVelocity[i] = 0;
				trgNb->FlowVelocity[idNb] = 0;
				continue;
			}
			/// 左右边界之间在固壁边界条件下不发生物质流动补偿
			if (vert1->IsOnLeftRightBoundary && vert2->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowVelocity[i] = 0;
				trgNb->FlowVelocity[idNb] = 0;
				continue;
			}
			///////////////////
			{
				// 与该条边对应的圆台面面积
				sideLength = CalcLength(Vec2DSub(trg->CycledPoses[i], trg->CycledPoses[(i + 1) % 3]));
				vertOpp = trg->Vertices[(i + 2) % 3];
				vertNbOpp = trgNb->Vertices[(idNb + 2) % 3];
				///////////////////////////////////////////////////////////////////////////////////// 第一步，计算三角形边上的流出速率
				///// 方案四
				moveDistance = (trg->FlowVelocity[i] + 0.5 * trg->FlowAcc[i] * halfDeltT) * halfDeltT;
				trg->FlowVelocity[i] += trg->FlowAcc[i] * halfDeltT;
				//double soundTime = 2.0 * trg->Area / sideLength / trg->SoundVelocity;
				//double soundTimeNb = 2.0 * trgNb->Area / sideLength / trgNb->SoundVelocity;
				//trg->FlowVelocity[i] -= trg->FlowVelocity[i] * max(0.0001, min(0.9999, (5.0 * halfDeltT / min(soundTime, soundTimeNb))));
				factorDecrease = max(trg->ViscCoeff / (trg->Area / 3), trgNb->ViscCoeff / (trgNb->Area / 3));
				trg->FlowVelocity[i] -= trg->FlowVelocity[i] * max(0.0001, min(0.9999, halfDeltT * factorDecrease));
				///// ..
				trgNb->FlowVelocity[idNb] = -trg->FlowVelocity[i];//让邻域三角形的物质流动速度同步
				//////////////////// 流动密度根据流动方向来取，如果是流出，就取为本三角形的密度。（是为了避免小密度三角形的流出量太大，造成负密度）
				if (moveDistance > 0) flowDensity = trg->Density;
				else flowDensity = trgNb->Density;
				//////////////////// 算出物质流动量
				massFlowOut = 0.5 * moveDistance * sideLength * flowDensity;
			}
			///////////////////////////////////////////////////////////////////////////////////// 第二步，物质流动引起的三角形内能和质量变化
			/////////////////////////////////////////////// 2.1 由物质流动携带内能引起的内能流动
			{
				if (massFlowOut >= 0)
				{
					internalEnergy = trg->InternalEnergy;
					trg->InternalEnergy -= (massFlowOut / trg->Mass) * internalEnergy;//注意，这里的内能修改要放在下面的质量修改语句之前。
					trgNb->InternalEnergy += (massFlowOut / trg->Mass) * internalEnergy;
				}
				else
				{
					internalEnergyNb = trgNb->InternalEnergy;
					trg->InternalEnergy -= (massFlowOut / trgNb->Mass) * internalEnergyNb;
					trgNb->InternalEnergy += (massFlowOut / trgNb->Mass) * internalEnergyNb;
				}
				/////////////////////////////////////////////// 2.2 物质流动造成的单元体积变化做功引起的内能流动。（20170507）
				{					
					//double energyFlowOut = 0.5 * ((trg->Pressure + trg->ViscosTensor.TraceAverage()) * massFlowOut / trg->Density
					//                              + (trgNb->Pressure + trgNb->ViscosTensor.TraceAverage()) * massFlowOut / trgNb->Density);
					energyFlowOut = 0.5 * (-(trg->Pressure + trg->Viscos) * trg->Area * log(1 - massFlowOut / trg->Mass)
						+ (trgNb->Pressure + trgNb->Viscos) * trgNb->Area * log(1 + massFlowOut / trgNb->Mass));
					//double energyFlowOut = 0.5 * (-(trg->Pressure + trg->ViscosTensor.TraceAverage() + trg->Stress.TraceAverage()) * trg->Volumn * log(1 - massFlowOut / trg->Mass)
					//                              + (trgNb->Pressure + trgNb->ViscosTensor.TraceAverage() + trgNb->Stress.TraceAverage()) * trgNb->Volumn * log(1 + massFlowOut / trgNb->Mass));					
					trg->InternalEnergy -= energyFlowOut;
					trgNb->InternalEnergy += energyFlowOut;
				}
			}
			/////////////////////////////////////////////// 2.3 修正动量和动能补偿（采用了新的动量流动算法 2020-8-10）
			{
				//////////////////// S1 记录原来动能和原来的动量
				initialKineticEnergy = 0.5 * vertOpp->Mass * CalcLengthSquar(vertOpp->Velocity)
										+ 0.5 * vertNbOpp->Mass * CalcLengthSquar(vertNbOpp->Velocity);
				initialMomentumOpp = CValueMultVec2D(vertOpp->Mass, vertOpp->Velocity);
				initialMomentumNbOpp = CValueMultVec2D(vertNbOpp->Mass, vertNbOpp->Velocity);
				
				//////////////////// S2 物质流动引起的格点质量变化
				trg->Mass -= massFlowOut;
				trgNb->Mass += massFlowOut;
				vertOpp->Mass -= massFlowOut / 3.0;
				vertNbOpp->Mass += massFlowOut / 3.0;
				
				// printf("%s, %d\n", __FUNCTION__, __LINE__);
				///////////// 2020-8-26
				if ((vert1->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall) || (vert1->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall) ||
					(vert2->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall) || (vert2->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall) ||
					(vertOpp->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall) || (vertOpp->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall) ||
					(vertNbOpp->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall) || (vertNbOpp->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall) ||
					MF_Method == 0)
				// if (MF_Method == 0)
				{
					// printf("%s, %d\n", __FUNCTION__, __LINE__);
					//////////////////// S3 根据动量守恒修正格点速度
					//// 中点的动量变化量
					momentumFlowOut = Vec2DMultCValue(Vec2DAdd(vert1->Velocity, vert2->Velocity), massFlowOut / 2.0);
					newMomentumOpp = Vec2DSub(initialMomentumOpp, Vec2DDivideCValue(momentumFlowOut, 3.0));
					newMomentumNbOpp = Vec2DAdd(initialMomentumNbOpp, Vec2DDivideCValue(momentumFlowOut, 3.0));
					vertOpp->Velocity = Vec2DDivideCValue(newMomentumOpp, vertOpp->Mass);
					vertNbOpp->Velocity = Vec2DDivideCValue(newMomentumNbOpp, vertNbOpp->Mass);	

					//////////////////// S4 根据动能守恒修正单元内能 
					deltMassVert = massFlowOut / 3.0;
					deltMomentumVert = Vec2DDivideCValue(momentumFlowOut, 3.0);
					deltKineticEnergy = 0.5 * (
						Dot(Vec2DSub(CValueMultVec2D(2, initialMomentumOpp), deltMomentumVert), deltMomentumVert) / (vertOpp->Mass - deltMassVert)
						- deltMassVert / (vertOpp->Mass - deltMassVert) * Dot(initialMomentumOpp, initialMomentumOpp) / vertOpp->Mass
						+ Dot(Vec2DSub(CValueMultVec2D(-2, initialMomentumNbOpp), deltMomentumVert), deltMomentumVert) / (vertNbOpp->Mass + deltMassVert)
						+ deltMassVert / (vertNbOpp->Mass + deltMassVert) * Dot(initialMomentumNbOpp, initialMomentumNbOpp) / vertNbOpp->Mass
						);
					trg->InternalEnergy += 0.5 * deltKineticEnergy;
					trgNb->InternalEnergy += 0.5 * deltKineticEnergy;
				}
				else
				{
					// printf("%s, %d\n", __FUNCTION__, __LINE__);
					//////////////////// S3 根据动量守恒修正格点速度
					// 参考论文 ” Matter flow method for alleviating checkerboard oscillations in triangular mesh SGH Lagrangian simulation“，Li Zhao, Bo Xiao, et al.
					lambda = Vec2DMultCValue(Vec2DAdd(vert1->Velocity, vert2->Velocity), massFlowOut / 6.0);
					a = 1.0 / vertOpp->Mass + 1.0 / vertNbOpp->Mass;
					b = CValueMultVec2D(2, Vec2DSub(CValueMultVec2D(1.0 / vertNbOpp->Mass, initialMomentumNbOpp), CValueMultVec2D(1.0 / vertOpp->Mass, initialMomentumOpp)));
					c = (1.0 / vertOpp->Mass - 1.0 / (vertOpp->Mass + massFlowOut / 3.0)) * CalcLengthSquar(initialMomentumOpp) 
						+ (1.0 / vertNbOpp->Mass - 1.0 / (vertNbOpp->Mass - massFlowOut / 3.0)) * CalcLengthSquar(initialMomentumNbOpp);
					// double rr = CalcLengthSquar(b) / (4*a*a) - c / a;

					if(MF_Method == 1){
						momentumFlowOut = CalculateMoveMomentumMatterFlow1(a, b, c, lambda);
					}
					if(MF_Method == 2){
						momentumFlowOut = CalculateMoveMomentumMatterFlow2(a, b, c);
					}

					newMomentumOpp = Vec2DSub(initialMomentumOpp, momentumFlowOut);
					newMomentumNbOpp = Vec2DAdd(initialMomentumNbOpp, momentumFlowOut);
					vertOpp->Velocity = Vec2DDivideCValue(newMomentumOpp, vertOpp->Mass);
					vertNbOpp->Velocity = Vec2DDivideCValue(newMomentumNbOpp, vertNbOpp->Mass);
				}
			}

			/////////////////////////////////////////////// 重设强度类型的因变量
			SetStrengthTypeDependentVariables(trg);
			SetStrengthTypeDependentVariables(trgNb);
			
		}
	}
}

#if 0
/// <summary>
/// “物质流动补偿”。
/// </summary>
/// <param name="halfDeltT"></param>
void MatterFlowEvolve(double deltT)
{
	int TrgsArrLen = MeshObj.TrgsArrLen;
	Triangle Trgs = MeshObj.Trgs;
	BoxBoundaryCondition TopBottomBoundaryCondition = MeshObj.TopBottomBoundaryCondition;
	/////////////////////////////////////////////////////////////////////////////////////////////
	
	CalculateMatterFlowAcc(); //计算边上的物质流动加速度（要放在计算完格点受力之后） 


	double halfDeltT = deltT * 0.5;//由于每条边要被遍历两次，因此将时间减半
	/////////////////////////////////////////////////////////////////////////////// 一边循环，一边就修改了单元和格点的质量、能量、速度
	int i, j, k;
	Triangle trg = NULL, trgNb = NULL;
	Vertex vert1, vert2, vertOpp, vertNbOpp;
	double massFlowOut;
	int idNb = -1;
	double sideLength, moveDistance, factorDecrease, flowDensity;
	double internalEnergy, internalEnergyNb, energyFlowOut, initialKineticEnergy, deltMassVert, deltKineticEnergy;
	vec2D initialMomentumOpp, initialMomentumNbOpp, momentumFlowOut, newMomentumOpp, newMomentumNbOpp, deltMomentumVert;
	
	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		if (trg->IsDead) continue;
		//////////////////////////////
		if (trg->MaterialId <= 0) continue;//真空和墙壁单元不计算物质流动
		//////////////////////////////////////////// 遍历三角形的每条边
		for (i = 0; i < 3; i++)
		{
			///////////////////////////////////////////////////////////////////////////////////// 找到该条边的相关邻域三角形、顶点和边信息
			trgNb = trg->NeighbourTrgs[i];			
			if(trgNb == NullTriangle) continue;
			
			/// 找到邻域三角形相关的顶点号
			for (idNb = -1, j = 0; j < 3; j++)
			{
				if (trgNb->NeighbourTrgs[j] == trg)
				{
					idNb = j;
					break;
				}
			}
			///////////////////
			if (trgNb->MaterialId != trg->MaterialId)// 不同种材料之间不发生物质流动补偿
			{
				trg->FlowVelocity[i] = 0;
				trgNb->FlowVelocity[idNb] = 0;
				continue;
			}
			/////////////////// 
			vert1 = trg->Vertices[i];
			vert2 = trg->Vertices[(i + 1) % 3];
			/// 上下底边之间在固壁边界条件下不发生物质流动补偿
			if (vert1->IsOnUpDownBoundary && vert2->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowVelocity[i] = 0;
				trgNb->FlowVelocity[idNb] = 0;
				continue;
			}
			/// 左右边界之间在固壁边界条件下不发生物质流动补偿
			if (vert1->IsOnLeftRightBoundary && vert2->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowVelocity[i] = 0;
				trgNb->FlowVelocity[idNb] = 0;
				continue;
			}
			///////////////////
			{
				// 与该条边对应的圆台面面积
				sideLength = CalcLength(Vec2DSub(trg->CycledPoses[i], trg->CycledPoses[(i + 1) % 3]));
				vertOpp = trg->Vertices[(i + 2) % 3];
				vertNbOpp = trgNb->Vertices[(idNb + 2) % 3];
				///////////////////////////////////////////////////////////////////////////////////// 第一步，计算三角形边上的流出速率
				///// 方案四
				moveDistance = (trg->FlowVelocity[i] + 0.5 * trg->FlowAcc[i] * halfDeltT) * halfDeltT;
				trg->FlowVelocity[i] += trg->FlowAcc[i] * halfDeltT;
				//double soundTime = 2.0 * trg->Area / sideLength / trg->SoundVelocity;
				//double soundTimeNb = 2.0 * trgNb->Area / sideLength / trgNb->SoundVelocity;
				//trg->FlowVelocity[i] -= trg->FlowVelocity[i] * max(0.0001, min(0.9999, (5.0 * halfDeltT / min(soundTime, soundTimeNb))));
				factorDecrease = max(trg->ViscCoeff / (trg->Area / 3), trgNb->ViscCoeff / (trgNb->Area / 3));
				trg->FlowVelocity[i] -= trg->FlowVelocity[i] * max(0.0001, min(0.9999, halfDeltT * factorDecrease));
				///// ..
				trgNb->FlowVelocity[idNb] = -trg->FlowVelocity[i];//让邻域三角形的物质流动速度同步
				//////////////////// 流动密度根据流动方向来取，如果是流出，就取为本三角形的密度。（是为了避免小密度三角形的流出量太大，造成负密度）
				if (moveDistance > 0) flowDensity = trg->Density;
				else flowDensity = trgNb->Density;
				//////////////////// 算出物质流动量
				massFlowOut = 0.5 * moveDistance * sideLength * flowDensity;
			}
			///////////////////////////////////////////////////////////////////////////////////// 第二步，物质流动引起的三角形内能和质量变化
			/////////////////////////////////////////////// 2.1 由物质流动携带内能引起的内能流动
			if (trg->MaterialId > 0)
			{
				if (massFlowOut >= 0)
				{
					internalEnergy = trg->InternalEnergy;
					trg->InternalEnergy -= (massFlowOut / trg->Mass) * internalEnergy;//注意，这里的内能修改要放在下面的质量修改语句之前。
					trgNb->InternalEnergy += (massFlowOut / trg->Mass) * internalEnergy;
				}
				else
				{
					internalEnergyNb = trgNb->InternalEnergy;
					trg->InternalEnergy -= (massFlowOut / trgNb->Mass) * internalEnergyNb;
					trgNb->InternalEnergy += (massFlowOut / trgNb->Mass) * internalEnergyNb;
				}
				/////////////////////////////////////////////// 2.2 物质流动造成的单元体积变化做功引起的内能流动。（20170507）
				{						
					//double energyFlowOut = 0.5 * ((trg->Pressure + trg->ViscosTensor.TraceAverage()) * massFlowOut / trg->Density
					//                              + (trgNb->Pressure + trgNb->ViscosTensor.TraceAverage()) * massFlowOut / trgNb->Density);
					energyFlowOut = 0.5 * (-(trg->Pressure + trg->Viscos) * trg->Area * log(1 - massFlowOut / trg->Mass)
						+ (trgNb->Pressure + trgNb->Viscos) * trgNb->Area * log(1 + massFlowOut / trgNb->Mass));
					//double energyFlowOut = 0.5 * (-(trg->Pressure + trg->ViscosTensor.TraceAverage() + trg->Stress.TraceAverage()) * trg->Volumn * log(1 - massFlowOut / trg->Mass)
					//                              + (trgNb->Pressure + trgNb->ViscosTensor.TraceAverage() + trgNb->Stress.TraceAverage()) * trgNb->Volumn * log(1 + massFlowOut / trgNb->Mass));
					
					trg->InternalEnergy -= energyFlowOut;
					trgNb->InternalEnergy += energyFlowOut;
				}
			}
			/////////////////////////////////////////////// 2.3 修正动量和动能补偿（采用了新的动量流动算法 20171031）
			if (trg->MaterialId > 0)
			{
				//////////////////// 记录原来动能和原来的动量
				initialKineticEnergy = 0.5 * vertOpp->Mass * CalcLengthSquar(vertOpp->Velocity)
										+ 0.5 * vertNbOpp->Mass * CalcLengthSquar(vertNbOpp->Velocity);
				initialMomentumOpp = CValueMultVec2D(vertOpp->Mass, vertOpp->Velocity);
				initialMomentumNbOpp = CValueMultVec2D(vertNbOpp->Mass, vertNbOpp->Velocity);
				//////////////////// 物质流动引起的格点质量变化
				trg->Mass -= massFlowOut;
				trgNb->Mass += massFlowOut;
				vertOpp->Mass -= massFlowOut / 3.0;
				vertNbOpp->Mass += massFlowOut / 3.0;
				{
					//////////////////// 根据动量守恒修正格点速度
					momentumFlowOut = Vec2DMultCValue(Vec2DAdd(vert1->Velocity, vert2->Velocity), massFlowOut / 2.0);
					newMomentumOpp = Vec2DSub(initialMomentumOpp, Vec2DDivideCValue(momentumFlowOut, 3.0));
					newMomentumNbOpp = Vec2DAdd(initialMomentumNbOpp, Vec2DDivideCValue(momentumFlowOut, 3.0));
					vertOpp->Velocity = Vec2DDivideCValue(newMomentumOpp, vertOpp->Mass);
					vertNbOpp->Velocity = Vec2DDivideCValue(newMomentumNbOpp, vertNbOpp->Mass);
					//////////////////// 根据动能守恒修正单元内能
					{
						deltMassVert = massFlowOut / 3.0;
						deltMomentumVert = Vec2DDivideCValue(momentumFlowOut, 3.0);
						deltKineticEnergy = 0.5 * (
							Dot(Vec2DSub(CValueMultVec2D(2, initialMomentumOpp), deltMomentumVert), deltMomentumVert) / (vertOpp->Mass - deltMassVert)
							- deltMassVert / (vertOpp->Mass - deltMassVert) * Dot(initialMomentumOpp, initialMomentumOpp) / vertOpp->Mass
							+ Dot(Vec2DSub(CValueMultVec2D(-2, initialMomentumNbOpp), deltMomentumVert), deltMomentumVert) / (vertNbOpp->Mass + deltMassVert)
							+ deltMassVert / (vertNbOpp->Mass + deltMassVert) * Dot(initialMomentumNbOpp, initialMomentumNbOpp) / vertNbOpp->Mass
							);
						trg->InternalEnergy += 0.5 * deltKineticEnergy;
						trgNb->InternalEnergy += 0.5 * deltKineticEnergy;
					}
				}
				//double newKineticEnergy = 0.5 * vertOpp->Mass * pow(vertOpp->Velocity.CalcLengthSquar()
				//                        + 0.5 * vertNbOpp->Mass * pow(vertNbOpp->Velocity.CalcLengthSquar();
				//trg->InternalEnergy += 0.5 * (initialKineticEnergy - newKineticEnergy);
				//trgNb->InternalEnergy += 0.5 * (initialKineticEnergy - newKineticEnergy);
			}
			else
			{
				//////////////////// 物质流动引起的格点质量变化
				trg->Mass -= massFlowOut;
				trgNb->Mass += massFlowOut;
				vertOpp->Mass -= massFlowOut / 3.0;
				vertNbOpp->Mass += massFlowOut / 3.0;
			}
			/////////////////////////////////////////////// 2.5 重设强度类型的因变量
			SetStrengthTypeDependentVariables(trg);
			SetStrengthTypeDependentVariables(trgNb);
		}
	}
}
#endif


/// <summary>
/// 设置强度类型的因变量，不包括周期坐标和面积。
/// </summary>
/// <param name="trg"></param>
void SetStrengthTypeDependentVariables(Triangle trg)
{
	int i;
	/////////////////////////////////////////////////////////////////////// 1. 更新密度、磁场应力张量
	trg->Density = trg->Mass / trg->Area;
	/////////////////////////////////////////////////////////////////////// 2. 根据状态方程计算压强和温度（20200628不需要二分法搜索温度了）
	if (trg->MaterialId > 0)
	{
		double specifiEnergyOfCell = trg->InternalEnergy / trg->Mass;
		EOSCalcPressureTemperature(trg->MaterialId, trg->Density, specifiEnergyOfCell, &trg->Pressure, &trg->Temperature);
	}

	/////////////////////////////////////////////////////////////////////// 计算三角形的最小高、声速
	double maxSideLength = 0.0;
	double youngSum = 0.0;
	/////////////////////////// 计算最小高
	if (trg->MaterialId >= 0)
	{
		double sideLength[3];
		for (i = 0; i < 3; i++)
		{
			sideLength[i] = CalcLength(Vec2DSub(trg->CycledPoses[i], trg->CycledPoses[(i + 1) % 3]));
		}
		maxSideLength = max(max(sideLength[0], sideLength[1]), sideLength[2]);
		{
			double * height = trg->Heights;;
			for (i = 0; i < 3; i++)
			{
				height[i] = 2 * trg->Area / sideLength[i];
			}
			trg->MinHeight = min(min(height[0], height[1]), height[2]);
		}
	}
	/////////////////////////// 计算总模量、弹性模量、声速备用
	if (trg->MaterialId >= 0)
	{
		if (trg->MaterialId > 0)
		{
			/////////////////////// 计算声速
			{
				/////// 第一步计算总模量
				/// 由流体声速计算体模量
				double soundVelocityFluid = CalcSoundVelocity(trg->MaterialId, trg->Density, trg->Temperature);
				double volumeYoung = soundVelocityFluid * soundVelocityFluid * trg->Density;
				/// 总模量
				youngSum = (volumeYoung);
			}
			/////// 第二步计算总模量确定的声速
			{
				double soundVelocity = sqrt(youngSum / trg->Density);
				trg->SoundVelocitySum = max(1e-20, soundVelocity);
			}
		}
		else
		{
			double soundVelocity = CalcSoundVelocity(trg->MaterialId, trg->Density, trg->Temperature);
			trg->SoundVelocitySum = max(1e-20, soundVelocity);
		}
	}	
	
	/////////////////////////////////////////////////////////////////////// 计算粘性系数和粘性应力。
	///////////////////////// 标量粘性	
	if (trg->MaterialId >= 0)
	{
		////////////// 计算应变率
		double strainRatio = CalcStrainRate(trg->CycledPoses[0], trg->CycledPoses[1], trg->CycledPoses[2],
			trg->Vertices[0]->Velocity, trg->Vertices[1]->Velocity, trg->Vertices[2]->Velocity);
		/// 计算物质流动对应变率张量的贡献
		if (trg->MaterialId > 0 && HAVE_MF)
		{
			// 计算物质流动造成的物质流出率
			double massFlowOutRatio = 0.0;
			for (i = 0; i < 3; i++)
			{
				double sideLength = CalcLength(Vec2DSub(trg->CycledPoses[i], trg->CycledPoses[(i + 1) % 3]));
				double flowDensity;
				{
					if (trg->FlowVelocity[i] >= 0) flowDensity = trg->Density;
					else{
						if(trg->NeighbourTrgs[i] != NullTriangle){
							flowDensity = trg->NeighbourTrgs[i]->Density;
						}else{
							printf("%s,%d, trg->NeighbourTrgs[%d] == NullTriangle\n", __FUNCTION__, __LINE__, i);	
						}
					}
				}
				massFlowOutRatio += 0.5 * trg->FlowVelocity[i] * sideLength * flowDensity;
			}
			// printf("物质流动对应变率张量的贡献 massFlowOutRatio = %e\n", massFlowOutRatio);
			// 物质流动对应变率张量的贡献
			strainRatio += (0.5 * massFlowOutRatio / trg->Mass);
		}
		////////////// 确定粘性系数（取动态人工粘性系数和物理粘性系数之间的较大值）
		{
			/// 采用 G. Scovazzi, JCP 231 (2012) 8029-8069 中的第 (50) 式。
			double artiViscCoeff = 2.0 * maxSideLength * maxSideLength * (-strainRatio);
			artiViscCoeff = max(artiViscCoeff, 1.0 * maxSideLength * trg->SoundVelocitySum);

			trg->ViscCoeff = max(artiViscCoeff, GetViscousCoefficient(trg->MaterialId));
			trg->ViscCoeff *= ViscCoeffTimes; // 2020.7.4
		}
		////////////// 从应变率和粘性系数计算粘性应力
		trg->Viscos = -trg->ViscCoeff * strainRatio * trg->Density;

		/////////////////////////////////////////////////////////////////////// 流体内应力之和。包括：压强、弹性应力、粘性应力。
		if (trg->MaterialId > 0)
		{
			trg->TotalFluidStress = CValueToTensor2D(trg->Pressure + trg->Viscos);//CValueAddTensor2D(trg->Pressure + trg->Viscos, trg->Stress);
		}
		else
		{
			trg->TotalFluidStress = CValueToTensor2D(trg->Pressure + trg->Viscos);//trg->Pressure * UnitTensor + trg->ViscosTensor.TraceAverage() * UnitTensor;
		}
	}
}

/// <summary>
/// 设置该三角形的所有因变量：面积，密度，压强,粘性力，标准长度.
/// </summary>
/// <param name="trg"></param>
void SetAllDependentVariables(Triangle trg)
{
	trg->CycledPoses[0] = trg->Vertices[0]->Pos;
	trg->CycledPoses[1] = trg->Vertices[1]->Pos;
	trg->CycledPoses[2] = trg->Vertices[2]->Pos;
	///////////////////////////////////// 计算面积
	trg->Area = CalcTrgArea(trg->CycledPoses[0], trg->CycledPoses[1], trg->CycledPoses[2]);
	// printf("trg->Area = %f\n", trg->Area);
	///////////////////////////////////// 设置强度类型的因变量
	SetStrengthTypeDependentVariables(trg);
}


/// <summary>
/// 设置所有三角形的所有因变量：面积，密度，压强,粘性力，标准长度.
/// </summary>
void SetAllDependentVariablesOfTrgs()
{
	int i, TrgsArrLen = MeshObj.TrgsArrLen;
	Triangle Trgs = MeshObj.Trgs,  trg;
	/////////////////////////////////////////////////////// 重设面积，密度，压强,粘性力，标准长度这些因变量
	for (i = 0; i < TrgsArrLen; i++)
	{
		trg = &Trgs[i];
		if (trg->IsDead) continue;
		//////////////////////////////
		SetAllDependentVariables(trg);
	}
}


/// <summary>
/// 从每个三角单元属性确定时间步长
/// </summary>
/// <returns></returns>
double DetermineDeltT()
{
	double deltTMinGlobal = 1e20;//预先给一个较大的初始值
	int deltTNameIDGlobal = -1;
	int i, TrgsArrLen = MeshObj.TrgsArrLen, t;
	Triangle Trgs = MeshObj.Trgs;
	//////////////////////////////////////////////////////////////////// 然后遍历三角形取最小的deltT
	int k;
	Triangle trg;
	double deltTMinOfThisTrg;
	int deltTNameIDOfThisTrg;
	double deltTSoundVelocity, deltTViscous, deltTVertexVelocity, deltTVertexAcc, deltTMatterFlow, deltTMatterFlowAcc;
	double maxHeightChangeRate, maxHeighChangeRateByAcc;
	double maxFlowRate, matterFlowRate[3];
	double maxFlowRateByAcc, matterFlowRateByAcc[3];
	///////////////////////////////////////

	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		deltTMinOfThisTrg = 1e20;//预先给一个较大的初始值
		deltTNameIDOfThisTrg = -1;
		//////////////////////////////
		if (trg->IsDead) continue;
		//////////////////////////////
		if (trg->MaterialId <= 0) continue;
		/////////////////////////////////////////////////////

		///////////////////////////////////////////////// (1) 计算声速决定的时间步长
		deltTSoundVelocity = timeStepSafeFactor * trg->MinHeight / trg->SoundVelocitySum;
		///////////////////////////////////////////////// (2) 计算粘性决定的时间步长, [G.Scovazzi2012JCP]
		deltTViscous = timeStepSafeFactor * (trg->MinHeight * trg->MinHeight / trg->ViscCoeff);
		///////////////////////////////////////////////// (3) 计算三角形顶点速度和加速度决定的时间步长
		{
			maxHeightChangeRate = CalcTrgMaxHeightChangeRate(
				Vec2DSub(trg->CycledPoses[1], trg->CycledPoses[0]),
				Vec2DSub(trg->CycledPoses[2], trg->CycledPoses[0]),
				Vec2DSub(trg->Vertices[1]->Velocity, trg->Vertices[0]->Velocity),
				Vec2DSub(trg->Vertices[2]->Velocity, trg->Vertices[0]->Velocity)
				);
			maxHeighChangeRateByAcc = CalcTrgMaxHeightChangeRate(
				Vec2DSub(trg->CycledPoses[1], trg->CycledPoses[0]),
				Vec2DSub(trg->CycledPoses[2], trg->CycledPoses[0]),
				Vec2DSub(Vec2DDivideCValue(trg->Vertices[1]->Force, trg->Vertices[1]->Mass), Vec2DDivideCValue(trg->Vertices[0]->Force, trg->Vertices[0]->Mass)),
				Vec2DSub(Vec2DDivideCValue(trg->Vertices[2]->Force, trg->Vertices[2]->Mass), Vec2DDivideCValue(trg->Vertices[0]->Force, trg->Vertices[0]->Mass))
				);
			/// 从最大高度变化率计算时间步长
			if (maxHeightChangeRate < 1e-20)
			{
				deltTVertexVelocity = deltTMinOfThisTrg;
			}
			else
			{
				deltTVertexVelocity = timeStepSafeFactor * (1.0 / maxHeightChangeRate);
			}
			/// 从加速度决定的最大高度变化率计算时间步长
			if (maxHeighChangeRateByAcc < 1e-20)
			{
				deltTVertexAcc = deltTMinOfThisTrg;
			}
			else
			{
				deltTVertexAcc = timeStepSafeFactor * sqrt(2.0 / maxHeighChangeRateByAcc);
			}
		}
		///////////////////////////////////////////////// (4) 计算物质流动速度决定的时间步长
		{
			{
				for (i = 0; i < 3; i++)
				{
					matterFlowRate[i] = fabs(trg->FlowVelocity[i] / trg->Heights[i]);
				}
				maxFlowRate = max(max(matterFlowRate[0], matterFlowRate[1]), matterFlowRate[2]);
			}
			/// 从流动速率计算时间步长
			if (maxFlowRate < 1e-20)
			{
				deltTMatterFlow = deltTMinOfThisTrg;
			}
			else
			{
				deltTMatterFlow = timeStepSafeFactor * (0.333 / maxFlowRate);
			}
		}
		///////////////////////////////////////////////// (5) 计算物质流动加速度决定的时间步长（20190410增）
		{
			{
				for (i = 0; i < 3; i++)
				{
					matterFlowRateByAcc[i] = sqrt(fabs(3.0 / 2 * trg->FlowAcc[i] / trg->Heights[i]));
				}
				maxFlowRateByAcc = max(max(matterFlowRateByAcc[0], matterFlowRateByAcc[1]), matterFlowRateByAcc[2]);
			}
			/// 从加速度引起的流动速率计算时间步长
			if (maxFlowRateByAcc < 1e-20)
			{
				deltTMatterFlowAcc = deltTMinOfThisTrg;
			}
			else
			{
				deltTMatterFlowAcc = timeStepSafeFactor / maxFlowRateByAcc;
				// deltTMatterFlowAcc = 0.1 * timeStepSafeFactor / maxFlowRateByAcc;
			}
		}

		///////////////////////////////////////////////// 从五种时间步长中取最小的
		if (deltTSoundVelocity < deltTMinOfThisTrg)
		{
			deltTMinOfThisTrg = deltTSoundVelocity;
			deltTNameIDOfThisTrg = 0;
		}
		if (deltTViscous < deltTMinOfThisTrg)
		{
			deltTMinOfThisTrg = deltTViscous;
			deltTNameIDOfThisTrg = 1;
		}
		if (deltTVertexAcc < deltTMinOfThisTrg)
		{
			deltTMinOfThisTrg = deltTVertexAcc;
			deltTNameIDOfThisTrg = 2;
		}
		if (deltTVertexVelocity < deltTMinOfThisTrg)
		{
			deltTMinOfThisTrg = deltTVertexVelocity;
			deltTNameIDOfThisTrg = 3;
		}
		if (deltTMatterFlow < deltTMinOfThisTrg)
		{
			deltTMinOfThisTrg = deltTMatterFlow;
			deltTNameIDOfThisTrg = 4;
		}
		if (deltTMatterFlowAcc < deltTMinOfThisTrg)
		{
			deltTMinOfThisTrg = deltTMatterFlowAcc;
			deltTNameIDOfThisTrg = 5;
		}

		//////////////////////////////////////////////
		if (deltTMinOfThisTrg < deltTMinGlobal)
		{
			deltTMinGlobal = deltTMinOfThisTrg;
			deltTNameIDGlobal = deltTNameIDOfThisTrg;
		}
	}
	////////////////////////////////////////////////////////////////////
	return (deltTMinGlobal);
}


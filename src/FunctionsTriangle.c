/*! \file  FunctionsTriangle.c
 *
 *  \brief Define some functions for Triangle structures 
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"


/// <summary>
/// 用于Split网格重分操作创建新三角形时的初始化
/// </summary>
void ClearTriangle(Triangle trg)
{
    trg->IsDead = false;
    trg->InternalEnergy = 0.0;
    trg->Mass = 0.0;
    trg->FlowAcc[0] = 0.0;
    trg->FlowAcc[1] = 0.0;
    trg->FlowAcc[2] = 0.0;
    trg->FlowVelocity[0] = 0.0;
    trg->FlowVelocity[1] = 0.0;
    trg->FlowVelocity[2] = 0.0;
    /////////////////////// 注意，千万不能把trg的Index清零。
}
/// <summary>
/// 用于网格重分时，将新三角形的主要内容（不是全部内容）复制给旧三角形
/// </summary>
/// <param name="trgSource"></param>
/// <param name="trgTarget"></param>
void CopyMainContentsTo(Triangle trgSource, Triangle trgTarget)
{
	int i;
    for(i=0; i<3; i++)
		trgTarget->Vertices[i] = trgSource->Vertices[i];
    for(i=0; i<3; i++)
		trgTarget->NeighbourTrgs[i] = trgSource->NeighbourTrgs[i];
    trgTarget->Mass = trgSource->Mass;
    trgTarget->MaterialId = trgSource->MaterialId;
    for(i=0; i<3; i++)
		trgTarget->FlowVelocity[i] = trgSource->FlowVelocity[i];
    if (trgSource->MaterialId > 0)
    {
        trgTarget->InternalEnergy = trgSource->InternalEnergy;
    }
}
/// <summary>
/// 设置三角形的三个顶点。
/// </summary>
void SetTrgVertices(Vertex verts[3], Vertex vert0, Vertex vert1, Vertex vert2)
{
	verts[0] = vert0;
	verts[1] = vert1;
	verts[2] = vert2;
}
/// <summary>
/// 设置三角形的三个顶点。
/// </summary>
void SetTrgVerticesByArray(Vertex verts[3], Vertex vertsSource[3])
{
	verts[0] = vertsSource[0];
	verts[1] = vertsSource[1];
	verts[2] = vertsSource[2];
}
/// <summary>
/// 设置三角形的三个邻域。
/// </summary>
void SetTrgNeighbours(Triangle trgNbs[3], Triangle trgNb0, Triangle trgNb1, Triangle trgNb2)
{
	trgNbs[0] = trgNb0;
	trgNbs[1] = trgNb1;
	trgNbs[2] = trgNb2;
}
/// <summary>
/// 设置三角形的三个顶点和三个邻域。
/// </summary>
void SetTrgVertsAndNbs(Triangle trg, Vertex vert1, Vertex vert2, Vertex vert3, Triangle trg1, Triangle trg2, Triangle trg3)
{
	Vertex * Vertices = trg->Vertices;
	Triangle * NeighbourTrgs = trg->NeighbourTrgs;
	////////////////////////
    Vertices[0] = vert1;
    Vertices[1] = vert2;
    Vertices[2] = vert3;
    NeighbourTrgs[0] = trg1;
    NeighbourTrgs[1] = trg2;
    NeighbourTrgs[2] = trg3;
}

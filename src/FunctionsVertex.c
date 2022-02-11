/*! \file  FunctionsVertex.c
 *
 *  \brief Define some functions for Vertex structures 
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"

/// When the entire full array of vertices is created, the index is initialized and does not need to be changed later
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

/// Clear vertex struct
void ClearVertex(Vertex vert)
{
    vert->IsDead = false;
    
	vert->IsOnUpDownBoundary 	= false;
    vert->IsOnUpBoundary 		= false;
    vert->IsOnDownBoundary 		= false;
	
	vert->IsOnLeftRightBoundary = false;
	vert->IsOnLeftBoundary		= false;
	vert->IsOnRightBoundary 	= false;
	
	vert->Mass = 0;
    ////////////////////////////////// Temporary auxiliary variables
	vert->Force = MakeVec2D(0,0);
    vert->ForceForMatterFlow = MakeVec2D(0,0);
}

bool VertexArrayContains(Vertex verts[], int vertsLen, Vertex ele)
{
	int i;
	for(i=0; i<vertsLen; i++)
	{
		if(verts[i] == ele)
		{
			return(true);
		}
	}
	return(false);
}


/// Set dead node
void SetDeadVertex(Vertex vert)
{
	vert->IsDead = true; // set dead flag
	MeshObj.DeadVertIdsArr[MeshObj.DeadVertIdsArrLen++] = vert->Index;
}

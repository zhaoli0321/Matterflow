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
    ////////////////////////////////// 临时辅助变量
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
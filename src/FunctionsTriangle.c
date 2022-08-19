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


/// When the entire full array of triangle is created, the index is initialized and does not need to be changed later
void CreateTrgsArray(struct triangle * * trgsArrOut)
{
	struct triangle * trgsArr = (struct triangle *)calloc(MAX_VERTS * 2, sizeof(struct triangle));
	int i;
	for (i = 0; i < MAX_VERTS * 2; i++)
	{
		trgsArr[i].Index = i;
	}
	(*trgsArrOut) = trgsArr;
}


/// Set the three vertices of the triangle.
void SetTrgVertices(Vertex verts[3], Vertex vert0, Vertex vert1, Vertex vert2)
{
	verts[0] = vert0;
	verts[1] = vert1;
	verts[2] = vert2;
}

/// Set the three vertices of the triangle.
void SetTrgVerticesByArray(Vertex verts[3], Vertex vertsSource[3])
{
	verts[0] = vertsSource[0];
	verts[1] = vertsSource[1];
	verts[2] = vertsSource[2];
}

/// Set the three neighbor triangles of the triangle.
void SetTrgNeighbours(Triangle trgNbs[3], Triangle trgNb0, Triangle trgNb1, Triangle trgNb2)
{
	trgNbs[0] = trgNb0;
	trgNbs[1] = trgNb1;
	trgNbs[2] = trgNb2;
}

/// Set the three vertices and three neighbors of the triangle.
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


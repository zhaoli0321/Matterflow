/*! \file  LoadMesh.c
 *
 *  \brief Define some load mesh functions
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"


void SetNodeFromFile(NodeFromFile node, double x, double y, double vX, double vY)
{
	node->Pos = MakeVec2D(x,y);
	node->Vel = MakeVec2D(vX,vY);
}

void SetTriangleFromFile(TriangleFromFile trg, int node1, int node2, int node3)
{
    trg->NodeIDs[0] = node1;
    trg->NodeIDs[1] = node2;
    trg->NodeIDs[2] = node3;
}

/// Read node data from the input file nodeListFile.
/// File format: Indices    Coordinate_x    Coordinate_y    Velocity_x   Velocity_y
void LoadNodes(char nodeListFile[], struct nodeFromFile nodesArray[], int * nodesArrLenOut)
{
	char line[1001], words[5][101];
    int nodesArrLen = 0;
    ////////////////////////////////////////////////////////// 1) open a file
    FILE * fr = MyOpenFile(nodeListFile, "r");
	////////////////////////////////////////////////////////// 2) skip first line
	if(SearchNextInFile(fr, "Velocity_y") == false)
	{
		ERR("Error reading first line of node file!\n");
	}
    ////////////////////////////////////////////////////////// 3) read data
    while (true)
    {
		fgets(line, 1000, fr);
		if(feof(fr))break;
        /////////////////////////////////////////
        words[4][0] = '\0';
		sscanf(line, "%s%s%s%s%s", words[0], words[1], words[2], words[3], words[4]);
        if (strcmp(words[4], "")!=0)
        {
            double x = StringToDouble(words[1]);
            double y = StringToDouble(words[2]);
            double vX = StringToDouble(words[3]);
            double vY = StringToDouble(words[4]);
            SetNodeFromFile(&nodesArray[nodesArrLen++], x, y, vX, vY);
        }
    }
    if (nodesArrLen == 0)
    {
        ERR("Error: The file does not contain valid node data!\n");
    }
    ////////////////////////////////////////////////////////// 4) close the file
    fclose(fr);
	(*nodesArrLenOut) = nodesArrLen;
}

/// Read triangular cell data from input file cellListFile
/// Compatible with or without magnetic field Bz (Tesla)
/// File format: Indices   Node1   Node2   Node3   MaterialID  Pressure    Temperature
void LoadCells(char cellListFile[], struct triangleFromFile trgsArray[], int * trgsArrLenOut)
{
	char line[1001], words[8][101];
    int trgsArrLen = 0;
    bool isHaveBz = false;
    ////////////////////////////////////////////////////////// 1) open a file
    FILE * fr = MyOpenFile(cellListFile, "r");
	////////////////////////////////////////////////////////// 2) skip first line
    // Compatible with or without magnetic field Bz (Tesla)
    fgets(line, 1000, fr);
    if (ContainStr(line, "Bz(Tesla)"))
	{
		isHaveBz = true;
	}else if(ContainStr(line, "Temperature"))
	{
		isHaveBz = false;
	}else{
        ERR("Error reading first line of cell file!\n");
    }

    ////////////////////////////////////////////////////////// 3) read data
    while (true)
    {
		fgets(line, 1000, fr);
		if(feof(fr))break;
        /////////////////////////////////////////
        words[6][0] = '\0';
        words[7][0] = '\0';
        if (isHaveBz)
        {
            sscanf(line, "%s%s%s%s%s%s%s%s", words[0], words[1], words[2], words[3], words[4], words[5], words[6], words[7]);
        }else{
            sscanf(line, "%s%s%s%s%s%s%s", words[0], words[1], words[2], words[3], words[4], words[5], words[6]);
        }
        if (isHaveBz)
        {
            if (strcmp(words[7], "")!=0)
            {
                int id1 = StringToInt(words[1]) - 1;
                int id2 = StringToInt(words[2]) - 1;
                int id3 = StringToInt(words[3]) - 1;
                TriangleFromFile trg = &trgsArray[trgsArrLen++] ;
                SetTriangleFromFile(trg, id1, id2, id3);
                ///////////
                {
                    int materialId = StringToInt(words[4]);
                    if (materialId > 0)
                    {
                        trg->MaterialId = ConvertUserInputMatNumToSystemMatID(materialId);
                    }
                    else
                    {
                        trg->MaterialId = materialId;
                    }
                }
                ///////////
                if (trg->MaterialId > 0)
                {
                    {
                        trg->Pressure    = StringToDouble(words[5]) * 1e-6; // Convert the unit 10^5Pa of the input file to the standard unit 100GPa in the program
                        trg->Temperature = StringToDouble(words[6]);
                    }
                }
            }
        }
        else
        {
            if (strcmp(words[6], "")!=0)
            {
                int id1 = StringToInt(words[1]) - 1;
                int id2 = StringToInt(words[2]) - 1;
                int id3 = StringToInt(words[3]) - 1;
                TriangleFromFile trg = &trgsArray[trgsArrLen++] ;
                SetTriangleFromFile(trg, id1, id2, id3);
                ///////////
                {
                    int materialId = StringToInt(words[4]);
                    if (materialId > 0)
                    {
                        trg->MaterialId = ConvertUserInputMatNumToSystemMatID(materialId);
                    }
                    else
                    {
                        trg->MaterialId = materialId;
                    }
                }
                ///////////
                if (trg->MaterialId > 0)
                {
                    {
                        trg->Pressure    = StringToDouble(words[5]) * 1e-6; // Convert the unit 10^5Pa of the input file to the standard unit 100GPa in the program
                        trg->Temperature = StringToDouble(words[6]);
                    }
                }
            }
        }
    }
    if (trgsArrLen == 0)
    {
        ERR("Error: The file does not contain valid cell data!\n");
    }
    ////////////////////////////////////////////////////////// 4) close the file
    fclose(fr);
    (*trgsArrLenOut) = trgsArrLen;
}

/// Build triangular neighborhood relations from a set of nodes and triangular cells.
void ConstructTrgNeighbours(struct nodeFromFile nodes[], int nodesLen, struct triangleFromFile trgs[], int trgsLen)
{
	int i, j, k, t;
    /// Step 1: Build the neighborhood triangle for each node
	/// initialization
	for (i = 0; i < nodesLen; i++)
	{
		nodes[i].TrgIDsLen = 0;
	}
    /// Traverse the triangular elements and build the neighborhood triangles of each node
    for (i = 0; i < trgsLen; i++)
    {
        int * nodeIds = trgs[i].NodeIDs;
        nodes[nodeIds[0]].TrgIDs[nodes[nodeIds[0]].TrgIDsLen++] = i;
        nodes[nodeIds[1]].TrgIDs[nodes[nodeIds[1]].TrgIDsLen++] = i;
        nodes[nodeIds[2]].TrgIDs[nodes[nodeIds[2]].TrgIDsLen++] = i;
    }
    
    /// Step 2: Judging from the common triangle of the two nodes of each edge as adjacent elements
    for (i = 0; i < trgsLen; i++)
    {
        TriangleFromFile trg = &trgs[i];
        for (j = 0; j < 3; j++)
        {
            NodeFromFile node1 = &nodes[trg->NodeIDs[j]];
            NodeFromFile node2 = &nodes[trg->NodeIDs[(j + 1) % 3]];
            int* trgIDs1 = node1->TrgIDs;
            int* trgIds2 = node2->TrgIDs;
            //////////////////////////////////////////////
            int neighbourTrgID = -1;
            for (k=0; k<node1->TrgIDsLen; k++)
            {
                int trgId1 = trgIDs1[k];
                if (trgId1 == i) continue;
                //////////////////////
                for (t=0; t<node2->TrgIDsLen; t++)
                {
                    int trgId2 = trgIds2[t];
                    if (trgId2 == trgId1)
                    {
                        neighbourTrgID = trgId2;
                        break;
                    }
                }
                if (neighbourTrgID != -1) break;
            }
            trg->NeighbourTrgIDs[j] = neighbourTrgID;
        }
        // printf("trg[%d]->NeighbourTrgIDs = {%d, %d, %d}\n", i, trg->NeighbourTrgIDs[0], trg->NeighbourTrgIDs[1],trg->NeighbourTrgIDs[2]);
    }
}

/// Generate standard grid points and elements required in the program from modeling file grid points and elements
void ConvertModelMeshToStandard(struct nodeFromFile nodesFromFile[], int totNodes, struct triangleFromFile trgsFromFile[], int totTrgs)
{
	int i,j;
	struct vertex * verts;
	struct triangle * trgs;
    /////////////////////////////////////////////////////////////
    MeshObj.VertsArrLen = totNodes;
    verts = MeshObj.Vertices;
    MeshObj.TrgsArrLen = totTrgs;
    trgs = MeshObj.Trgs;
    /////////////////////////////////////////////////////////////
    for (i = 0; i < totNodes; i++)
    {
        Vertex vert = &verts[i];
        ClearVertex(vert);
        vert->Pos = nodesFromFile[i].Pos;
        vert->Velocity = nodesFromFile[i].Vel;       
    }
    /// Set three nodes of a triangle
    for (i = 0; i < totTrgs; i++)
    {
        int* nodeIDs = trgsFromFile[i].NodeIDs;
        ////////////////////////////////////////
        trgs[i].Vertices[0] = &verts[nodeIDs[0]];
        trgs[i].Vertices[1] = &verts[nodeIDs[1]];
        trgs[i].Vertices[2] = &verts[nodeIDs[2]];
        ////////////////////////////////////////       
    }

    ///////////////////////////////////////////////////////////// Set the physical properties of the triangle
    int* neighbourIDs;
    for (i = 0; i < totTrgs; i++)
    {
        Triangle trg = &trgs[i];

        /// Create adjacent triangle elements
        neighbourIDs = trgsFromFile[i].NeighbourTrgIDs;
        for (j = 0; j < 3; j++)
        {
            if (neighbourIDs[j] != -1)
            {
                trgs[i].NeighbourTrgs[j] = &trgs[neighbourIDs[j]];
            }
            else
            {
                trgs[i].NeighbourTrgs[j] = NullTriangle;
            }
        }

        ////////////////////////////////////////////// Calculate the area of a triangular element
		double area;
        area = CalcTrgArea(trg->Vertices[0]->Pos, trg->Vertices[1]->Pos, trg->Vertices[2]->Pos);
		
        ////////////////////////////////////////////// Set basic physical quantities such as density, temperature, etc.
        if (trgsFromFile[i].MaterialId > 0)
        {
            ////////////////////////////// 
            Material material = &MatParasList[trgsFromFile[i].MaterialId];
            {
                trg->MaterialId = trgsFromFile[i].MaterialId;
                //////////////////////////////////////////////
                trg->Temperature = trgsFromFile[i].Temperature;
                trg->Density = EOSCalcDensity(trg->MaterialId, trgsFromFile[i].Pressure, trgsFromFile[i].Temperature);
            }
            ////////////////////////////// 
            trg->Mass = trg->Density * area;
            {
                double pressure, specificEnergy;
                EOSCalcPressureEnergy(trg->MaterialId, trg->Density, trg->Temperature, &pressure, &specificEnergy);
                trg->InternalEnergy = trg->Mass * specificEnergy;
            }
        }
        else if (trgsFromFile[i].MaterialId == 0) // If it is a vacuum unit, it is described by a linear elastic fluid with minimal density
        {
            trg->MaterialId = 0;
            trg->Density = MatParasList[0].NormalDensity;
            trg->Mass = trg->Density * area;
        }
        else // In the case of the Wall Material element, the linear elastic fluid description is also used.
        {
            trg->MaterialId = trgsFromFile[i].MaterialId;
            trg->Density = MatWall[-trg->MaterialId].NormalDensity;
            trg->Mass = trg->Density * area;
		}
    }
    ///////////////////////////////////////////////////////////// Set mass of the node
    for (i = 0; i < totTrgs; i++)
    {
        Triangle trg = &trgs[i];
        //////////////////////////////
		{
			double oneThirdTrgMass = 1 / 3.0 * trg->Mass;
			int t;
			for (t = 0; t < 3; t++ )
			{
				Vertex vertex = trg->Vertices[t];
				vertex->Mass += oneThirdTrgMass;
			}
		}
    }
}

/// Create the properties of the grid boundary 
void SetUpDownLeftRightBoundaryVertices(struct vertex vertices[], int vertsArrLen)
{
	int i;
    //////////////////////////////////// First get minX, minY, tinyLength , so that it is convenient to judge whether the node is on the boundary.
    double minX, minY, maxX, maxY, tinyLength;
    {
        minX = maxX = vertices[0].Pos.X;
        minY = maxY = vertices[0].Pos.Y;
        for (i = 0; i < vertsArrLen; i++ )
        {
            Vertex vert = &vertices[i];
            //////////////////////////
			{
				double posX = vert->Pos.X;
				if (minX > posX) minX = posX;
				else if (maxX < posX) maxX = posX;
			}
            /////
			{
				double posY = vert->Pos.Y;
				if (minY > posY) minY = posY;
				else if (maxY < posY) maxY = posY;
			}
        }
        tinyLength = 0.0001 * sqrt((maxX - minX) * (maxY - minY) / vertsArrLen);
    }
    //////////////////////////////////// judge whether the node is on the boundary.
    for (i = 0; i < vertsArrLen; i++)
    {
        Vertex vert = &vertices[i];
        //////////////////////////
        if (fabs(vert->Pos.X - minX) < tinyLength) vert->IsOnLeftRightBoundary = true;       
        /////
        if (fabs(vert->Pos.Y - minY) < tinyLength) vert->IsOnUpDownBoundary = true;


        if (fabs(vert->Pos.X - maxX) < tinyLength) vert->IsOnLeftRightBoundary = true;
        if (fabs(vert->Pos.Y - maxY) < tinyLength) vert->IsOnUpDownBoundary = true;

        if (fabs(vert->Pos.X - minX) < tinyLength) vert->IsOnLeftBoundary = true;  
        if (fabs(vert->Pos.X - maxX) < tinyLength) vert->IsOnRightBoundary = true;

        if (fabs(vert->Pos.Y - minY) < tinyLength) vert->IsOnDownBoundary = true;
        if (fabs(vert->Pos.Y - maxY) < tinyLength) vert->IsOnUpBoundary = true;   
    }
}

/// Shift the coordinates - the bottom left corner is zeroed. And determine the mesh width and height.
void TranslatePositions(struct nodeFromFile nodes[], int totNodes, struct triangleFromFile trgs[], int totTrgs, double * widthOut, double * heightOut)
{
	int i;
    /// Determine the translation amount according to the maximum and minimum coordinate values
    double minX, maxX, minY, maxY;
    minX = maxX = nodes[0].Pos.X;
    minY = maxY = nodes[0].Pos.Y;
    for(i= 0; i<totNodes; i++)
    {
        NodeFromFile node = &nodes[i];
        /////
		{
			double posX = node->Pos.X;
			if (minX > posX) minX = posX;
			else if (maxX < posX) maxX = posX;
		}
        /////
		{
			double posY = node->Pos.Y;
			if (minY > posY) minY = posY;
			else if (maxY < posY) maxY = posY;
		}
    }
    //////////////////////// Coordinate translation (the bottom left corner is zeroed)
	{
		vec2D translation = MakeVec2D(-minX, -minY);
		for ( i = 0; i < totNodes; i++)
		{
			NodeFromFile node = &nodes[i];
			Vec2DSelfAdd(&node->Pos, translation);
		}
	}
    //////////////////////// Domain width and height
    (*widthOut) = maxX - minX;
    (*heightOut) = maxY - minY;
}

/// Read mesh file and build mesh
void SetMeshFromMeshFiles(char nodesFileName[], char cellsFileName[])
{
    // 1) Import node and cell data
    struct nodeFromFile * nodesFromFile = (struct nodeFromFile *) malloc(sizeof(struct nodeFromFile) * MAX_VERTS);
    struct triangleFromFile * trgsFromFile = (struct triangleFromFile *) malloc(sizeof(struct triangleFromFile) * MAX_VERTS * 2);
    int totNodes, totTrgs;
    LoadNodes(nodesFileName, nodesFromFile, &totNodes);
    LoadCells(cellsFileName, trgsFromFile, &totTrgs);
    // 2) Shift the coordinates - the bottom left corner is zeroed. And determine the mesh width and height.
	TranslatePositions(nodesFromFile, totNodes, trgsFromFile, totTrgs, &MeshObj.Width, &MeshObj.Height);
    // 3) Build triangle beighborhood relations
    ConstructTrgNeighbours(nodesFromFile, totNodes, trgsFromFile, totTrgs);
    // 4) Generate standard grid points and elements required in the program from modeling files
    ConvertModelMeshToStandard(nodesFromFile, totNodes, trgsFromFile, totTrgs);
    // 5) Create the properties of the grid boundary 
    SetUpDownLeftRightBoundaryVertices(MeshObj.Vertices, MeshObj.VertsArrLen);
    /////////////////////////////////////////////////
	free(nodesFromFile);
	free(trgsFromFile);
}
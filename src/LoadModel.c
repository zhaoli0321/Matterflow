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


/// Read grid filenames, including grid point filenames and cell filenames
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
            ERR("Error reading LoadMeshFileNames: The file format does not meet the requirements!\n");
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

/// Read Computational Domain Boundary Conditions
void LoadBoxBoundaryCondition(FILE * fr, BoxBoundaryCondition * topBottomCondition, BoxBoundaryCondition * leftRightCondition)
{
	char line[1001], words[3][101];
	////////////////////////////////////////////////////////////////////
	(*topBottomCondition) = BoxBoundaryCondition_Wall;
	(*leftRightCondition) = BoxBoundaryCondition_Wall;
	////////////////////////////////////////////////////////////////////
	if(SearchNextInFile(fr, ">> Boundary Conditions of Computation Domain") == false)
	{
		ERR("Does not find \">> Boundary Conditions of Computation Domain.\"");
	}
	////////////////////////////////////////////////////////////////////
    while (true)
    {
        fgets(line, 1000, fr);
		if(feof(fr))
        {
            ERR("Error reading \"Boundary Conditions of Computation Domain\": The file format does not meet the requirements!\n");
        }
        /////////////////////////////////////////
        words[1][0] = '\0';
		sscanf(line, "%s%s%s", words[0], words[1], words[2]);
		if(strcmp(words[1], "=")==0)
		{
			if (strcmp(words[0], "TopBottomBoundaryCondition") == 0)
			{
				if (strcmp(words[2], "SolidWall") == 0)
				{
					(*topBottomCondition) = BoxBoundaryCondition_Wall;
				}
				else if (strcmp(words[2], "Period") == 0)
				{
					(*topBottomCondition) = BoxBoundaryCondition_Cycle;
					ERR("Periodic boundary conditions are currently not implemented!!!\n");
				}
				else
				{
					ERR("Error reading \"Boundary Conditions of Computation Domain\": The file format does not meet the requirements!\n");
				}
			}
			else if (strcmp(words[0], "LeftRightBoundaryCondition") == 0)
			{
				if (strcmp(words[2], "SolidWall") == 0)
				{
					(*leftRightCondition) = BoxBoundaryCondition_Wall;
				}
				else if (strcmp(words[2], "Period") == 0)
				{
					(*leftRightCondition) = BoxBoundaryCondition_Cycle;
					ERR("Periodic boundary conditions are currently not implemented!!!\n");
				}
				else
				{
					ERR("Error reading \"Boundary Conditions of Computation Domain\": The file format does not meet the requirements!\n");
				}
				break;
			}
		}
    }
}

/// Read Gravitational Acceleration
double ReadGravityAcceleration(FILE * fr)
{
	char line[1001], words[3][101];
	////////////////////////////////////////////////////////////////////
    double gravityAcc = 0.0;
	////////////////////////////////////////////////////////////////////
	if(SearchNextInFile(fr, ">> Gravitational Acceleration") == false)
	{
		ERR("Does not find \">> Gravitational Acceleration\"");
	}
	////////////////////////////////////////////////////////////////////
    while (true)
    {
        fgets(line, 1000, fr);
		if(feof(fr))
        {
			ERR("Error reading \"Gravitational Acceleration\": The file format does not meet the requirements!\n");
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

/// Read Time Interval of Data Output and Program End Time
void LoadTimeInterval(FILE * fr, double * timeInterval, double * timeEnd)
{
	char line[1001], words[3][101];
	////////////////////////////////////////////////////////////////////
    (*timeInterval) = 0.0;
    (*timeEnd) = 0.0;
	////////////////////////////////////////////////////////////////////
	if(SearchNextInFile(fr, ">> Time Interval of Data Output and Program End Time") == false)
	{
		ERR("Does not find \">> Time Interval of Data Output and Program End Time\"");
	}
	////////////////////////////////////////////////////////////////////
    while (true)
    {
		if(feof(fr))
        {
            ERR("Error reading \"Time Interval of Data Output and Program End Time\": The file format does not meet the requirements!\n");
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

/// Read the computational model from the input file modelFile.
void MeshObjectFromLoadModel(char * modelFile)
{
	FILE * fr;
    char nodesFileName[101];
    char cellsFileName[101];
    // 1) Open a file
    fr = MyOpenFile(modelFile, "r");
    // 2) Read the material list
    LoadMaterialListFromFile(fr);
    // 3) Read mesh
    LoadMeshFileNames(fr, nodesFileName, cellsFileName);
    SetMeshFromMeshFiles(nodesFileName, cellsFileName);
    // 4) Read calculation area boundary conditions, gravitational acceleration
    LoadBoxBoundaryCondition(fr, &MeshObj.TopBottomBoundaryCondition, &MeshObj.LeftRightBoundaryCondition);
    MeshObj.GravityFactor = ReadGravityAcceleration(fr);
    // 5) Set all dependent variables
    SetAllDependentVariablesOfTrgs();
    // 6) Set the modulus of the cavity according to the maximum sound velocity at the initial moment
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
    // 7) Time interval for reading data output
    LoadTimeInterval(fr, &TimeInterval, &TimeEnd);
    // 8) Close the file
    fclose(fr);
}
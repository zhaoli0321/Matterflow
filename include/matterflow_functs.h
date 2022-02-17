/*! \file  matterflow_functs.h
 *
 *  \brief Function decoration for the Matterflow package
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021-2022 by the CAEP-XTU team. All rights reserved.             
 *---------------------------------------------------------------------------------
 *
 *  \warning DO NOT EDIT!!! This file is automatically generated!
 */ 

#include "matterflow.h" 

/*-------- In file: FileAndStringReadWrite.c --------*/

FILE * MyOpenFile(char * fileName, char * format);

bool SearchNextInFile(FILE * fp, char * searchStr);

bool ContainStr(char * strSource, char * searchStr);

double StringToDouble(char * str);

int StringToInt(char * str);


/*-------- In file: FunctionsMaterial.c --------*/

void MakeVacuumMaterial(Material matVacuum);

void LoadMaterialListFromFile(FILE * fr);

int ConvertUserInputMatNumToSystemMatID(int matNum);

void EOSCalcPressureEnergy(int matId, double density, double temperature, double * pressure, double * specificEnergy);

void EOSCalcPressureTemperature(int matId, double density, double specificEnergy, double * pressure, double * temperature);

double EOSCalcDensity(int matId, double pressure, double temperature);

double CalcViscosStress(int matId, double strainRatio, double density, double cViscDynArti);

double GetViscousCoefficient(int matId);

double CalcSoundVelocity(int matId, double density, double temperature);

double CalcSpecificEnergyToTemperature(int matId, double density, double temperature);

tensor2D CalcViscosStressTensor(int matId, tensor2D strainRatioTensor, double density, double cViscDynArti);


/*-------- In file: FunctionsTriangle.c --------*/

void CreateTrgsArray(struct triangle * * trgsArrOut);

void SetTrgVertices(Vertex verts[3], Vertex vert0, Vertex vert1, Vertex vert2);

void SetTrgVerticesByArray(Vertex verts[3], Vertex vertsSource[3]);

void SetTrgNeighbours(Triangle trgNbs[3], Triangle trgNb0, Triangle trgNb1, Triangle trgNb2);

void SetTrgVertsAndNbs(Triangle trg, Vertex vert1, Vertex vert2, Vertex vert3, Triangle trg1, Triangle trg2, Triangle trg3);


/*-------- In file: FunctionsVertex.c --------*/

void CreateVerticesArray(struct vertex * * vertsArrOut);

void ClearVertex(Vertex vert);

bool VertexArrayContains(Vertex verts[], int vertsLen, Vertex ele);


/*-------- In file: GlobalParameters.c --------*/


/*-------- In file: LoadMesh.c --------*/

void SetNodeFromFile(NodeFromFile node, double x, double y, double vX, double vY);

void SetTriangleFromFile(TriangleFromFile trg, int node1, int node2, int node3);

void LoadNodes(char nodeListFile[], struct nodeFromFile nodesArray[], int * nodesArrLenOut);

void LoadCells(char cellListFile[], struct triangleFromFile trgsArray[], int * trgsArrLenOut);

void ConstructTrgNeighbours(struct nodeFromFile nodes[], int nodesLen, struct triangleFromFile trgs[], int trgsLen);

void ConvertModelMeshToStandard(struct nodeFromFile nodesFromFile[], int totNodes, struct triangleFromFile trgsFromFile[], int totTrgs);

void SetUpDownLeftRightBoundaryVertices(struct vertex vertices[], int vertsArrLen);

void TranslatePositions(struct nodeFromFile nodes[], int totNodes, struct triangleFromFile trgs[], int totTrgs, double * widthOut, double * heightOut);

void SetMeshFromMeshFiles(char nodesFileName[], char cellsFileName[]);


/*-------- In file: LoadModel.c --------*/

void LoadMeshFileNames(FILE * fr, char nodesFileName[], char cellsFileName[]);

void LoadBoxBoundaryCondition(FILE * fr, BoxBoundaryCondition * topBottomCondition, BoxBoundaryCondition * leftRightCondition);

double ReadGravityAcceleration(FILE * fr);

void LoadTimeInterval(FILE * fr, double * timeInterval, double * timeEnd);

void MeshObjectFromLoadModel(char * modelFile);


/*-------- In file: MeshObject.c --------*/

void ERR(char * prnt);

double GetWallTime(clock_t start_time, clock_t end_time);

double func(vec2D x, double a, vec2D b, double c);

vec2D CalculateMoveMomentumMatterFlow1(double a, vec2D b, double c, vec2D lambda);

vec2D CalculateMoveMomentumMatterFlow2(double a, vec2D b, double c);

void SetGravity();

void CalculateVertexForce();

void VelocityPosEdgeCondition();

void ForceEdgeCondition();

void DynamicEvolve(double deltT);

void CalculateMatterFlowAcc();

void MatterFlowEvolve(double deltT);

void SetStrengthTypeDependentVariables(Triangle trg);

void SetAllDependentVariables(Triangle trg);

void SetAllDependentVariablesOfTrgs();

double DetermineDeltT();


/*-------- In file: PostProcessing.c --------*/

void CreateDirectory(char *dirName);

void PrintNSInformation();

void PrintInfo();

void PrintMeshScale();

void WriteTecPlot();

double PressureAnalyticalSolution_TaylorGreen(double xc, double yc);

void VelocityAnalyticalSolution_TaylorGreen(double x, double y, double *vx, double *vy);

void TaylorGreenL2Norm();

void VolumeWeightedL1Norm();

void TaylorGreenOutPutBottomEdgePhysics();


/*-------- In file: Tensor2DMath.c --------*/

tensor2D Tensor2DAdd(tensor2D t1, tensor2D t2);

void Tensor2DSelfAdd(tensor2D * t1, tensor2D t2);

tensor2D Tensor2DSub(tensor2D t1, tensor2D t2);

tensor2D CValueAddTensor2D(double c, tensor2D t2);

tensor2D Tensor2DAddCValue(tensor2D t2, double c);

void Tensor2DSelfAddCValue(tensor2D * t1, double c);

tensor2D CValueMultTensor2D(double factor, tensor2D tensor);

tensor2D Tensor2DMultCValue(tensor2D tensor, double factor);

void Tensor2DSelfMultCValue(tensor2D * tensor, double factor);

vec2D Tensor2DMultVec(tensor2D tensor, vec2D vec);

tensor2D Tensor2DMult(tensor2D t1, tensor2D t2);

void SelfMultiply(tensor2D * t, double factor);

tensor2D Tensor2DDivideCValue(tensor2D tensor, double factor);

tensor2D DisplaceTensor(vec2D x2, vec2D x3, vec2D x2p, vec2D x3p);

double Trace(tensor2D t);

double TraceAverage(tensor2D t);

double GetTensor2DMaxValue(tensor2D t);

tensor2D CValueToTensor2D(double c);

void TensorDepartion(tensor2D t, tensor2D * Strain, tensor2D * Rotate);

tensor2D CalcStrainRateTensor(vec2D x1, vec2D x2, vec2D v1, vec2D v2);

double CalcStrainRate(vec2D x1, vec2D x2, vec2D x3, vec2D v1, vec2D v2, vec2D v3);

tensor2D Diagonalize(tensor2D t);


/*-------- In file: Vec2DMath.c --------*/

vec2D MakeVec2D(double r, double z);

void SelfRotate90(vec2D * v);

vec2D Rotate90(vec2D v);

vec2D Rotate270(vec2D v);

vec2D Vec2DAdd(vec2D v1, vec2D v2);

vec2D Vec2DAdd3(vec2D v1, vec2D v2, vec2D v3);

void Vec2DSelfAdd(vec2D * v1, vec2D v2);

vec2D Vec2DSub(vec2D v1, vec2D v2);

void Vec2DSelfSub(vec2D * v1, vec2D v2);

double Dist(vec2D pos1, vec2D pos2);

double CalcLength(vec2D vec);

double CalcLengthSquar(vec2D vec);

vec2D Vec2DMultCValue(vec2D vec, double ratio);

void Vec2DSelfMultCValue(vec2D * v1, double c);

vec2D CValueMultVec2D(double ratio, vec2D vec);

double Cross(vec2D v1, vec2D v2);

double Dot(vec2D v1, vec2D v2);

double CalcTrgArea(vec2D pos0, vec2D pos1, vec2D pos2);

vec2D Vec2DDivideCValue(vec2D vec, double ratio);

vec2D DirectVec(vec2D v1, vec2D v2);

vec2D UnitVec(vec2D vec);

double CalcTrgMaxHeightChangeRate(vec2D r1, vec2D r2, vec2D v1, vec2D v2);

/* End of matterflow_functs.h */

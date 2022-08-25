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

// Error function
void ERR(char * prnt)
{
    printf(prnt);
    getchar();
    exit(0);
}

// Get wall time
double GetWallTime(clock_t start_time, clock_t end_time)
{
	return ( (double) (end_time - start_time) / CLOCKS_PER_SEC );
}


// The Matter Flow method
/*
	min{ |lambda - p| }
	a<p, p> + <b, p> + c = 0
 */
vec2D CalculateMoveMomentumMatterFlow(double a, vec2D b, double c, vec2D lambda)
{
	vec2D x_O;			// x_O is the coordinates of the center of the circle
	double r, d, ratio; // r is the radius of the circle, d is the distance from the point lambda to x_O
	
	if(a < 0){
		printf("error a = %e < 0\n", a);
	}
	x_O = CValueMultVec2D(- 1.0 / (2*a), b);
	d = CalcLength( Vec2DSub(lambda, x_O) );
	double rr = CalcLengthSquar(b) / (4*a*a) - c / a;
	if(rr <= 0)
	{
		return (x_O);
	}else{
		r = sqrt(rr);
	} 

	if(d > 0)
	{
		ratio = r / d;
	}else{
		vec2D rn = {r, 0};
		return ( Vec2DAdd(x_O, rn) );
	}
	/////////////
	return ( Vec2DAdd(CValueMultVec2D(ratio, lambda), CValueMultVec2D(1-ratio, x_O)) );
}


/// Set gravity for each Vertex
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
		//////////////////////////////////////
		vert->Force = MakeVec2D(0, -MeshObj.GravityFactor * vert->Mass);
		vert->ForceForMatterFlow = ZeroVec2D;
	}
}


/// Calculate Nodal Force
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
		//////////////////////////////
		if (trg->MaterialId < 0) continue;
		///////////////////////////////////////////////////////////////////////////////////////////////// If it is a cell of matter, include the magnetic pressure
		{
			trgVertices = trg->Vertices;
			if (trg->MaterialId > 0)
			{
				////////////////////////////////////////////////////////////////// Calculation of Forces for Material Flow
				{						
					stress = CValueToTensor2D(trg->Pressure + trg->Viscos);					
					///////////////////////////////////// Divide the stress into three nodes
					for (i = 0; i < 3; i++)
					{
						/////////////////////////////// 
						side = Rotate270(Vec2DSub(trg->CycledPoses[(i + 1) % 3], trg->CycledPoses[i]));
						f = Tensor2DMultVec(stress, side);
						Vec2DSelfAdd(&trgVertices[i]->ForceForMatterFlow, Vec2DMultCValue(f, 0.5));
						Vec2DSelfAdd(&trgVertices[(i + 1) % 3]->ForceForMatterFlow, Vec2DMultCValue(f, 0.5));
					}
				}
				////////////////////////////////////////////////////////////////// Calculate full force
				{
					stress = trg->TotalFluidStress;
					///////////////////////////////////// Divide the stress into three nodes
					for (i = 0; i < 3; i++)
					{
						/////////////////////////////// 
						side = Rotate270(Vec2DSub(trg->CycledPoses[(i + 1) % 3], trg->CycledPoses[i]));
						/////////////////////////////// 
						f = Tensor2DMultVec(stress, side);
						Vec2DSelfAdd(&trgVertices[i]->Force, Vec2DMultCValue(f, 0.5));
						Vec2DSelfAdd(&trgVertices[(i + 1) % 3]->Force, Vec2DMultCValue(f, 0.5));
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////// If it is a vacuum cell
			else
			{
				stress = trg->TotalFluidStress;
				///////////////////////////////////// Divide the stress into three nodes
				for (i = 0; i < 3; i++)
				{
					/////////////////////////////// 
					side = Rotate270(Vec2DSub(trg->CycledPoses[(i + 1) % 3], trg->CycledPoses[i]));
					/////////////////////////////// 
					f = Tensor2DMultVec(stress, side);
					Vec2DSelfAdd(&trgVertices[i]->Force, Vec2DMultCValue(f, 0.5));
					Vec2DSelfAdd(&trgVertices[(i + 1) % 3]->Force, Vec2DMultCValue(f, 0.5));
				}
			}
		}
	}
}


/// Position/Velocity Boundary Conditions
void VelocityPosEdgeCondition()
{
	int VertsArrLen = MeshObj.VertsArrLen, TrgsArrLen = MeshObj.TrgsArrLen;
	Vertex Vertices = MeshObj.Vertices;
	Triangle Trgs = MeshObj.Trgs;
	BoxBoundaryCondition TopBottomBoundaryCondition = MeshObj.TopBottomBoundaryCondition;
	double Width  = MeshObj.Width;
	double Height = MeshObj.Height;
	/////////////////////////////////////////////////////////////////////////
	int k;
	Vertex vert;
	///////////////////////////////////////
	for (k = 0; k < VertsArrLen; k++)
	{
		vert = &Vertices[k];
		/////////// If the up/bottom edge is a solid wall boundary condition
		if (vert->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall)
		{
			vert->Velocity.Y = 0;
			if (vert->IsOnUpBoundary)	vert->Pos.Y = Height;	
			if (vert->IsOnDownBoundary)	vert->Pos.Y = 0.0;							
		}
		/////////// If the left/right edge is a solid wall boundary condition
		if (vert->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall)
		{
			vert->Velocity.X = 0;
			if (vert->IsOnRightBoundary)	vert->Pos.X = Width;				
			if (vert->IsOnLeftBoundary)		vert->Pos.X = 0.0;							
		}
	}

	////////////// Velocity conditions satisfied by points on wall elements
	int t;
	Triangle trg;
	///////////////////////////////////////
	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
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

/// Forced Boundary Conditions
void ForceEdgeCondition()
{
	int VertsArrLen = MeshObj.VertsArrLen, TrgsArrLen = MeshObj.TrgsArrLen;
	Vertex Vertices = MeshObj.Vertices;
	Triangle Trgs = MeshObj.Trgs;
	BoxBoundaryCondition TopBottomBoundaryCondition = MeshObj.TopBottomBoundaryCondition;
	int k;
	Vertex vert;
	///////////////////////////////////////

	for (k = 0; k < VertsArrLen; k++)
	{
		vert = &Vertices[k];
		///////////// The longitudinal force of the bottom node is zero
		if (vert->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall)
		{
			vert->Force.Y = 0;
			vert->ForceForMatterFlow.Y = 0;
		}
		/////////// The lateral force of the left and right boundary nodes is zero
		if (vert->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall)
		{
			vert->Force.X = 0;
			vert->ForceForMatterFlow.X = 0;
		}
	}

	///////////////////////////////////////
	int t;
	Triangle trg;
	///////////////////////////////////////
	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
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

/// Update the physical quantities at the next moment
void DynamicEvolve(double deltT)
{
	int VertsArrLen = MeshObj.VertsArrLen, TrgsArrLen = MeshObj.TrgsArrLen;
	Vertex Vertices = MeshObj.Vertices;
	Triangle Trgs = MeshObj.Trgs;
	/////////////////////////////////////////////////////// Update coordinates and velocity
	int k;
	Vertex vert;
	vec2D deltVelocity;
	///////////////////////////////////////
	for (k = 0; k < VertsArrLen; k++)
	{
		vert = &Vertices[k];
		//////////////////////////////////////
		{
			deltVelocity = Vec2DMultCValue(vert->Force, 1 / vert->Mass * deltT);
			// Record the displacement increment, which is used to calculate the internal energy increment below
			vert->DeltPos = Vec2DMultCValue(Vec2DAdd(vert->Velocity, CValueMultVec2D(0.5, deltVelocity)), deltT);
			Vec2DSelfAdd(&vert->Pos, vert->DeltPos);
			Vec2DSelfAdd(&vert->Velocity, deltVelocity);
		}	
	}

	/////////////////////////////////////////////////////// Update fluid internal energy (magnetic field stress does not contribute to fluid internal energy)
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
		//////////////////////////////
		if (trg->MaterialId > 0)
		{
			/// Convert the internal stress of this triangle to the force of the three nodes
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
				/////////////////////////////// ith side of the triangle
				side = Rotate270(Vec2DSub(trg->CycledPoses[(i + 1) % 3], trg->CycledPoses[i]));
				/////////////////////////////// Divide the stress into three nodes
				f = Tensor2DMultVec(stress, side);
				Vec2DSelfAdd(&fLats[i], Vec2DMultCValue(f, 0.5));
				Vec2DSelfAdd(&fLats[(i + 1) % 3], Vec2DMultCValue(f, 0.5));
			}
			/// Calculate the internal energy increment (stress work)
			{
				verts = trg->Vertices;
				trg->InternalEnergy += (-Dot(fLats[0], verts[0]->DeltPos) - Dot(fLats[1], verts[1]->DeltPos) - Dot(fLats[2], verts[2]->DeltPos));
			}
			/// Calculate the internal energy increment (source term) for the Taylor-Green vortex problem
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

/// Artificial heat flux: update the internal energy of the fluid
#if ThermalDiffusion
void DynamicEvolveThermalDiffusion(double deltT)
{
	int VertsArrLen = MeshObj.VertsArrLen, TrgsArrLen = MeshObj.TrgsArrLen;
	Vertex Vertices = MeshObj.Vertices;
	Triangle Trgs = MeshObj.Trgs;
	double Eps = 1E-20;

	int i, k;
	Triangle trg;
	double FluxRate;

	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		// The vacuum element is not allowed to have internal energy diffusion
		if (trg->MaterialId <= 0) continue; 
#if SALTZMAN			
		if (trg->MaterialId == 1) continue; // material == 1 是活塞单元
#endif	

		for (i = 0; i < 3; i++)
		{
			Triangle trgRight = trg->NeighbourTrgs[i];
			// If there is no element on the ith side of the triangle, skip
			if(trgRight == NullTriangle) continue;
#if SALTZMAN			
			if (trgRight->MaterialId == 1) continue; // material == 1 是活塞单元
#endif	
#if HalfMesh_SALTZMAN
			// 由于 Saltzman 的网格是镜像翻倍网格，
			// 为了降低计算量，程序中存的镜像翻倍的网格，但计算的时候只用一半的网格来计算
			if(trgRight->Index >= TrgsArrLen) continue; // 不流向另一半网格
#endif						
			/// The two end positions of the edge and the edge vector
			vec2D posEnd1 = trg->CycledPoses[i];
			vec2D posEnd2 = trg->CycledPoses[(i + 1) % 3];
			vec2D sideVec = Vec2DSub(posEnd2, posEnd1);
			double sideLength = CalcLength(sideVec);
			/// Distances from left and right center points to the edge
			double distLeft, distRight;
			/// Left distance
			{
				vec2D posCLeft = Vec2DDivideCValue(Vec2DAdd3(trg->CycledPoses[0], trg->CycledPoses[1], trg->CycledPoses[2]), 3.0);
				distLeft = fabs(Cross(Vec2DSub(posCLeft, posEnd1), sideVec)) / sideLength;
			}
			/// Right distance
			{
				vec2D posCRight = Vec2DDivideCValue(Vec2DAdd3(trgRight->CycledPoses[0], trgRight->CycledPoses[1], trgRight->CycledPoses[2]), 3.0);
				distRight = fabs(Cross(Vec2DSub(posCRight, posEnd1), sideVec)) / sideLength;
			}
			
			
			/// Internal energy densities in left and right triangles
			trg->InternalEnergyDensity = trg->InternalEnergy / trg->Mass; // Internal energy density in left triangle
			trgRight->InternalEnergyDensity = trgRight->InternalEnergy / trgRight->Mass; // Internal energy density in right triangle

			/// The thermal diffusion coefficient is equal to 0 and no heat flux occurs.
			if (fabs(trg->ThermalDiffusionCoeff) < Eps || fabs(trgRight->ThermalDiffusionCoeff) < Eps) FluxRate = 0.0;
			else 
			{
				double conductivityDistLeft = distLeft / trg->ThermalDiffusionCoeff; 
				double conductivityDistRight = distRight / trgRight->ThermalDiffusionCoeff;
				//printf("conductivityDistLeft = %e, conductivityDistRight = %e\n", conductivityDistLeft, conductivityDistRight);
				
				/// heat flux
				FluxRate = (trgRight->InternalEnergyDensity - trg->InternalEnergyDensity) * sideLength / (conductivityDistLeft + conductivityDistRight);
				if (isNAN(FluxRate)) FluxRate = 0.0;
			}

			// update the internal energy of element
			trg->InternalEnergy += deltT * FluxRate;
		}
	}
}
#endif

/// Calculate the acceleration of matter flow
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
		//////////////////////////////
		if (trg->MaterialId <= 0) continue; // Vacuum and wall cells do not use matter flow method
		//////////////////////////////////////////// 
		for (i = 0; i < 3; i++)
		{
			/////////////////////////////////////////////////////////////////////////////////////
			vert1 = trg->Vertices[i];
			vert2 = trg->Vertices[(i + 1) % 3];
			/// No matter flow compensation occurs between the upper and lower bases under the solid wall boundary condition
			if (vert1->IsOnUpDownBoundary && vert2->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowAcc[i] = 0.0;
				continue;
			}
			/// No matter flow compensation occurs between the left and right boundaries under the solid wall boundary condition
			if (vert1->IsOnLeftRightBoundary && vert2->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowAcc[i] = 0.0;
				continue;
			}
			///////////////////////////////////////////////////////////////////////////////////// Find the relevant neighbor triangles, vertices, and edge information for the edge
			trgNb = trg->NeighbourTrgs[i];			
			if(trgNb == NullTriangle) continue;
			
			/// Find the vertex number associated with the neighborhood triangle
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
			// No matter flow compensation occurs between different materials
			if (trgNb->MaterialId != trg->MaterialId) 
			{
				trg->FlowAcc[i] = 0.0;
				continue;
			}
			{
				///////////////////////////////////////////////////////////////////////////////////// Calculate edge-related information
				side = Rotate90(Vec2DSub(trg->CycledPoses[i], trg->CycledPoses[(i + 1) % 3]));
				sn = UnitVec(side);
				///////////////////////////////////////////////////////////////////////////////////// Calculate outflow rate acceleration on triangle sides
				//////////////////// Acceleration of two vertices and their average acceleration
				a1 = Vec2DDivideCValue(vert1->ForceForMatterFlow, vert1->Mass);
				a2 = Vec2DDivideCValue(vert2->ForceForMatterFlow, vert2->Mass);
				a12 = Dot(CValueMultVec2D(0.5, Vec2DAdd(a1, a2)), sn); 
				//////////////////// The acceleration of the assumed midpoint				
				force 	= CValueMultVec2D((trg->Pressure + trg->Viscos) * 0.5, side);
				forceNb = CValueMultVec2D((trgNb->Pressure + trgNb->Viscos) * (-0.5), side);					
				ac = Dot(Vec2DAdd(force, forceNb), sn) / ((trg->Mass + trgNb->Mass) / 3.0);
				//////////////////// 
				trg->FlowAcc[i] = ac - a12;
			}
		}
	}
}


/// Matter Flow Evolve
void MatterFlowEvolve(double deltT)
{	
	int TrgsArrLen = MeshObj.TrgsArrLen;
	Triangle Trgs = MeshObj.Trgs;
	BoxBoundaryCondition TopBottomBoundaryCondition = MeshObj.TopBottomBoundaryCondition;
	/////////////////////////////////////////////////////////////////////////////////////////////
	double halfDeltT = deltT * 0.5; // Cuts the time in half since each edge is traversed twice

	CalculateMatterFlowAcc();

	/////////////////////////////////////////////////////////////////////////////// Modified mass, energy, velocity of elements and nodes
	int i, j, k;
	Triangle trg = NULL, trgNb = NULL;
	Vertex vert1, vert2, vertOpp, vertNbOpp;
	double massFlowOut;
	int idNb = -1;
	double sideLength, moveDistance, C_diss, flowDensity;
	double internalEnergy, internalEnergyNb, energyFlowOut, initialKineticEnergy, deltMassVert, deltKineticEnergy;
	vec2D initialMomentumOpp, initialMomentumNbOpp, momentumFlowOut, newMomentumOpp, newMomentumNbOpp, deltMomentumVert;
	double a, c;
	vec2D b, lambda;

	for (k = 0; k < TrgsArrLen; k++)
	{
		trg = &Trgs[k];
		//////////////////////////////
		if (trg->MaterialId <= 0) continue; // Vacuum and wall cells do not use matter flow method
		//////////////////////////////////////////// 
		for (i = 0; i < 3; i++)
		{
			///////////////////////////////////////////////////////////////////////////////////// Find the relevant neighbor triangles, vertices, and edge information for the edge
			trgNb = trg->NeighbourTrgs[i];			
			if(trgNb == NullTriangle) continue;

			/// Find the vertex number associated with the neighborhood triangle
			for (idNb = -1, j = 0; j < 3; j++)
			{
				if (trgNb->NeighbourTrgs[j] == trg)
				{
					idNb = j;
					break;
				}
			}
			///////////////////
			// No matter flow compensation occurs between different materials
			if (trgNb->MaterialId != trg->MaterialId) 
			{
				trg->FlowVelocity[i] = 0;
				trgNb->FlowVelocity[idNb] = 0;
				continue;
			}
			/////////////////// 
			vert1 = trg->Vertices[i];
			vert2 = trg->Vertices[(i + 1) % 3];
			/// No matter flow compensation occurs between the upper and lower bases under the solid wall boundary condition
			if (vert1->IsOnUpDownBoundary && vert2->IsOnUpDownBoundary && MeshObj.TopBottomBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowVelocity[i] = 0;
				trgNb->FlowVelocity[idNb] = 0;
				continue;
			}
			/// No matter flow compensation occurs between the left and right boundaries under the solid wall boundary condition
			if (vert1->IsOnLeftRightBoundary && vert2->IsOnLeftRightBoundary && MeshObj.LeftRightBoundaryCondition == BoxBoundaryCondition_Wall)
			{
				trg->FlowVelocity[i] = 0;
				trgNb->FlowVelocity[idNb] = 0;
				continue;
			}
			///////////////////////////////////// Step 1: Calculate the outflow rate on the sides of the triangle
			{
				sideLength = CalcLength(Vec2DSub(trg->CycledPoses[i], trg->CycledPoses[(i + 1) % 3]));
				vertOpp = trg->Vertices[(i + 2) % 3];
				vertNbOpp = trgNb->Vertices[(idNb + 2) % 3];
				moveDistance = trg->FlowVelocity[i] * halfDeltT + 0.5 * trg->FlowAcc[i] * halfDeltT * halfDeltT;
				trg->FlowVelocity[i] += trg->FlowAcc[i] * halfDeltT;
				
				// An artificial dissipation factor C_diss is introduced in the evolution equation of the matter flow velocity.
				C_diss = (trg->ViscCoeff * trg->Density / trg->Area + trgNb->ViscCoeff * trgNb->Density / trgNb->Area) * sideLength * sideLength / ((trg->Mass + trgNb->Mass) / 3.0);
				trg->FlowVelocity[i] -= trg->FlowVelocity[i] * halfDeltT * C_diss;
				
				///
				trgNb->FlowVelocity[idNb] = -trg->FlowVelocity[i]; // Synchronize matter flow velocities in neighborhood triangles
				// The flow density is taken according to the flow direction. If it is outflow, it is taken as the density of this triangle.
				//  (It is to avoid too much outflow of small density triangles, resulting in negative density)
				if (moveDistance > 0) flowDensity = trg->Density;
				else flowDensity = trgNb->Density;
				//////////////////// Calculate matter outflow
				massFlowOut = 0.5 * moveDistance * sideLength * flowDensity;
			}
			///////////////////////////////////// Step 2: Changes in the internal energy and mass of the triangle caused by the flow of matter
			///////////////////////////////////// Step 2.1 Internal energy flow caused by internal energy carried by matter flow
			{
				if (massFlowOut >= 0)
				{
					// Note that the energy modification here should be placed before the quality modification statement below.
					internalEnergy = trg->InternalEnergy;
					trg->InternalEnergy -= (massFlowOut / trg->Mass) * internalEnergy; 
					trgNb->InternalEnergy += (massFlowOut / trg->Mass) * internalEnergy;
				}
				else
				{
					internalEnergyNb = trgNb->InternalEnergy;
					trg->InternalEnergy -= (massFlowOut / trgNb->Mass) * internalEnergyNb;
					trgNb->InternalEnergy += (massFlowOut / trgNb->Mass) * internalEnergyNb;
				}
				///////////////////////////////////// Step 2.2 The flow of internal energy due to the work done by the volume change of the unit caused by the flow of matter.
				{					
					energyFlowOut = 0.5 * ( + (trg->Pressure   + trg->Viscos  ) * massFlowOut / trg->Density
											+ (trgNb->Pressure + trgNb->Viscos) * massFlowOut / trgNb->Density);
					trg->InternalEnergy -= energyFlowOut;
					trgNb->InternalEnergy += energyFlowOut;
				}
			}
			///////////////////////////////////// Step 2.3 Corrected momentum and kinetic energy compensation
			{
				//////////////////// S1 Record old kinetic energy and momentum
				initialKineticEnergy = 0.5 * vertOpp->Mass * CalcLengthSquar(vertOpp->Velocity)
									 + 0.5 * vertNbOpp->Mass * CalcLengthSquar(vertNbOpp->Velocity);
				initialMomentumOpp = CValueMultVec2D(vertOpp->Mass, vertOpp->Velocity);
				initialMomentumNbOpp = CValueMultVec2D(vertNbOpp->Mass, vertNbOpp->Velocity);
				
				//////////////////// S2 Nodal mass change due to matter flow
				trg->Mass -= massFlowOut;
				trgNb->Mass += massFlowOut;
				vertOpp->Mass -= massFlowOut / 3.0;
				vertNbOpp->Mass += massFlowOut / 3.0;

				// printf("%s, %d\n", __FUNCTION__, __LINE__);
				//////////////////// S3 Correction of grid velocities
				lambda = Vec2DMultCValue(Vec2DAdd(vert1->Velocity, vert2->Velocity), massFlowOut / 6.0);
				a = 1.0 / vertOpp->Mass + 1.0 / vertNbOpp->Mass;
				b = CValueMultVec2D(2, Vec2DSub(CValueMultVec2D(1.0 / vertNbOpp->Mass, initialMomentumNbOpp), CValueMultVec2D(1.0 / vertOpp->Mass, initialMomentumOpp)));
				c = (1.0 / vertOpp->Mass - 1.0 / (vertOpp->Mass + massFlowOut / 3.0)) * CalcLengthSquar(initialMomentumOpp) 
					+ (1.0 / vertNbOpp->Mass - 1.0 / (vertNbOpp->Mass - massFlowOut / 3.0)) * CalcLengthSquar(initialMomentumNbOpp);

				// obtaining outflow momentum by solving optimization problem
				momentumFlowOut = CalculateMoveMomentumMatterFlow(a, b, c, lambda); 

				newMomentumOpp   = Vec2DSub(initialMomentumOpp, momentumFlowOut);
				newMomentumNbOpp = Vec2DAdd(initialMomentumNbOpp, momentumFlowOut);
				vertOpp->Velocity   = Vec2DDivideCValue(newMomentumOpp, vertOpp->Mass);
				vertNbOpp->Velocity = Vec2DDivideCValue(newMomentumNbOpp, vertNbOpp->Mass);
			}

			/////////////////////////////////////////////// Reset the dependent variables of the intensity type
			SetStrengthTypeDependentVariables(trg);
			SetStrengthTypeDependentVariables(trgNb);
		}
	}
}


/// Sets the dependent variables of type intensity, excluding periodic coordinates and area.
void SetStrengthTypeDependentVariables(Triangle trg)
{
	int i;
	/////////////////////////////////////////////////////////////////////// Update density, magnetic field stress tensors
	trg->Density = trg->Mass / trg->Area;
	/////////////////////////////////////////////////////////////////////// Calculate pressure and temperature from equation of state
	if (trg->MaterialId > 0)
	{
		double specifiEnergyOfCell = trg->InternalEnergy / trg->Mass;
		EOSCalcPressureTemperature(trg->MaterialId, trg->Density, specifiEnergyOfCell, &trg->Pressure, &trg->Temperature);
	}

	/////////////////////////////////////////////////////////////////////// Calculate the minimum height and speed of sound of a triangle
	double maxSideLength = 0.0;
	double youngSum = 0.0;
	/////////////////////////// Calculate the minimum height
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
	/////////////////////////// Calculate total modulus, elastic modulus, sound velocity
	if (trg->MaterialId >= 0)
	{
		if (trg->MaterialId > 0)
		{
			/////////////////////// Calculate the speed of sound
			{
				/////// S1: Calculate total modulus
				/// Calculate Bulk Modulus from Fluid Sound Velocity
				double soundVelocityFluid = CalcSoundVelocity(trg->MaterialId, trg->Density, trg->Temperature);
				double volumeYoung = soundVelocityFluid * soundVelocityFluid * trg->Density;
				/// total modulus
				youngSum = (volumeYoung);
			}
			/////// S2: Calculate the speed of sound determined by the total modulus
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
	
	/////////////////////////////////////////////////////////////////////// Calculate the viscous coefficient and viscous stress.
	///////////////////////// Scalar viscosity	
	if (trg->MaterialId >= 0)
	{
		////////////// Calculate strain rate
		double strainRatio = CalcStrainRate(trg->CycledPoses[0], trg->CycledPoses[1], trg->CycledPoses[2],
			trg->Vertices[0]->Velocity, trg->Vertices[1]->Velocity, trg->Vertices[2]->Velocity);
		
		/// Calculate the contribution of matter flow to the strain rate
		if (trg->MaterialId > 0 && HAVE_MF)
		{
			// Calculate the rate of matter outflow 
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
			// Contribution of matter flow to the strain rate
			strainRatio += massFlowOutRatio / trg->Mass;
		}
		
		/// Determine the viscosity coefficient
		double L = sqrt(trg->Area0); // Using constant length
		trg->ViscCoeff = ViscCoeffTimes * max(-strainRatio * L * L, trg->SoundVelocitySum * L);
		if (isNAN(trg->ViscCoeff)) trg->ViscCoeff = 0;

		/// Calculate viscous stress from strain rate and viscosity coefficient
		trg->Viscos = -trg->ViscCoeff * strainRatio * trg->Density;

		/// The sum of the internal stresses in the fluid. Including: pressure and viscous stress.
		trg->TotalFluidStress = CValueToTensor2D(trg->Pressure + trg->Viscos);

#if ThermalDiffusion
		double c1 = 0.5, c2 = 1.0;	
		double cs = trg->SoundVelocitySum; // speed of sound
		
		// The heat diffusion coefficient
		double Kappa = TDCoeffTimes * (c1 * L * L * fabs(strainRatio) + c2 * cs * L);
		trg->ThermalDiffusionCoeff = Kappa * trg->Density;
#endif
	}
	else
	{ 
		// No thermal diffusion in vacuum elements
		trg->ThermalDiffusionCoeff = 0; 
	}
}

/// Set all dependent variables of this triangle: area, density, pressure, viscous force.
void SetAllDependentVariables(Triangle trg)
{
	trg->CycledPoses[0] = trg->Vertices[0]->Pos;
	trg->CycledPoses[1] = trg->Vertices[1]->Pos;
	trg->CycledPoses[2] = trg->Vertices[2]->Pos;
	///////////////////////////////////// area
	trg->Area = CalcTrgArea(trg->CycledPoses[0], trg->CycledPoses[1], trg->CycledPoses[2]);
	// printf("trg->Area = %f\n", trg->Area);
	///////////////////////////////////// Set the dependent variables of type intensity,
	SetStrengthTypeDependentVariables(trg);
}


/// Set all dependent variables for all triangles: area, density, pressure, viscous force.
void SetAllDependentVariablesOfTrgs()
{
	int i, TrgsArrLen = MeshObj.TrgsArrLen;
	Triangle Trgs = MeshObj.Trgs,  trg;
	/////////////////////////////////////////////////////// Set all dependent variables
	for (i = 0; i < TrgsArrLen; i++)
	{
		trg = &Trgs[i];
		//////////////////////////////
		SetAllDependentVariables(trg);
	}
}


/// Determining the time-step size from each triangular element property
double DetermineDeltT()
{
	double deltTMinGlobal = 1e20; // give a large initial value
	int deltTNameIDGlobal = -1;
	int i, TrgsArrLen = MeshObj.TrgsArrLen, t;
	Triangle Trgs = MeshObj.Trgs;
	//////////////////////////////////////////////////////////////////// 
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
		deltTMinOfThisTrg = 1e20;
		deltTNameIDOfThisTrg = -1;
		//////////////////////////////
		if (trg->MaterialId <= 0) continue;
		///

		/// (1) Calculate the time-step determined by the speed of sound
		deltTSoundVelocity = timeStepSafeFactor * trg->MinHeight / trg->SoundVelocitySum;
		
		/// (2) Calculate the time-step determined by the viscosity
		if (trg->ViscCoeff == 0) deltTViscous = 1e20;
		else deltTViscous = timeStepSafeFactor * (trg->MinHeight * trg->MinHeight / (2*trg->ViscCoeff));
		
		/// (3-4) Calculate the time-step determined by the velocities and accelerations of the triangle vertices
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
			/// Calculate time-step from maximum rate of height change
			if (maxHeightChangeRate < 1e-20)
			{
				deltTVertexVelocity = deltTMinOfThisTrg;
			}
			else
			{
				deltTVertexVelocity = timeStepSafeFactor * (1.0 / maxHeightChangeRate);
			}
			/// Calculate the time step from the maximum rate of change of height determined by acceleration
			if (maxHeighChangeRateByAcc < 1e-20)
			{
				deltTVertexAcc = deltTMinOfThisTrg;
			}
			else
			{
				deltTVertexAcc = timeStepSafeFactor * sqrt(2.0 / maxHeighChangeRateByAcc);
			}
		}
		/// (5) Calculate the time-step determined by the velocity of the matter flow
		{
			{
				for (i = 0; i < 3; i++)
				{
					matterFlowRate[i] = fabs(trg->FlowVelocity[i] / trg->Heights[i]);
				}
				maxFlowRate = max(max(matterFlowRate[0], matterFlowRate[1]), matterFlowRate[2]);
			}
			/// 
			if (maxFlowRate < 1e-20)
			{
				deltTMatterFlow = deltTMinOfThisTrg;
			}
			else
			{
				deltTMatterFlow = timeStepSafeFactor * (0.333 / maxFlowRate);
			}
		}
		///////////////////////////////////////////////// (6) Calculate the time-step determined by the acceleration of the matter flow
		{
			{
				for (i = 0; i < 3; i++)
				{
					matterFlowRateByAcc[i] = sqrt(fabs(3.0 / 2 * trg->FlowAcc[i] / trg->Heights[i]));
				}
				maxFlowRateByAcc = max(max(matterFlowRateByAcc[0], matterFlowRateByAcc[1]), matterFlowRateByAcc[2]);
			}
			/// 
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

		///////////////////////////////////////////////// Take the smallest time-step
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
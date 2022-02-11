/*! \file  FunctionsMaterial.c
 *
 *  \brief Define some functions for Material structures 
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"



/// Create a vacuum material
void MakeVacuumMaterial(Material matVacuum)
{
	/////////////////////////////////////////////
	matVacuum->MaterialNumber = 0;
	strcpy(matVacuum->MaterialName, "Vacuum");
	///////////////////////
	matVacuum->EOSType = EOSTypeEnum_Linear;
	matVacuum->NormalDensity = 1e-8;    // Room temperature/pressure density 
	matVacuum->Kai = 1e-10;             // Bulk modulus
	matVacuum->SpecificHeat = 1.0;
}


/// Import material data list from file
void LoadMaterialListFromFile(FILE * fr)
{
	char line[1001], words[3][1001];
	////////////////////////////////////////////// 0th element put vacuum material
	MakeVacuumMaterial(&MatParasList[MatParasListLen++]);
	////////////////////////////////////////////// 
	if (SearchNextInFile(fr, ">> Material Settings") == false)
	{
		ERR("Does not find \">> Material Settings\"\n");
	}
	////////////////////////////////////////////// Read the material list
	{
		Material lastMatParas = NULL;
		int materialNumber = 1;
		bool relateMaterialIsNeed = false; // Flag if associated material needs to be created
		while (true)
		{
			///////////////////////////////////////// read a material
			fgets(line, 1000, fr);
			if (feof(fr))
			{
				ERR("Material format does not meet requirements\n");
			}
			if (ContainStr(line, "--------------------"))
			{
				////////////////////////////////////////// Create a new material
				MatParasListLen++;
				lastMatParas = &MatParasList[MatParasListLen - 1];
				lastMatParas->MaterialNumber = materialNumber++;
				relateMaterialIsNeed = false;
			}
			else if (ContainStr(line, "===================="))
			{
				break;
			}
			else
			{
				/////////////////////////////////////////
				words[1][0] = '\0';
				sscanf(line, "%s%s%s", words[0], words[1], words[2]);
				if (strcmp(words[1], "=") == 0)
				{
					if (strcmp(words[0], "MaterialName") == 0)
					{
						strcpy(lastMatParas->MaterialName, words[2]);
					}
					else if (strcmp(words[0], "EOSType") == 0)
					{
						int typeNum = StringToInt(words[2]);
						switch (typeNum)
						{
						case 1:
							lastMatParas->EOSType = EOSTypeEnum_Linear;
							break;
						case 2:
							lastMatParas->EOSType = EOSTypeEnum_IdealGas;
							break;
						}
					}
					else if (strcmp(words[0], "NormalDensity") == 0)
					{
						lastMatParas->NormalDensity = StringToDouble(words[2]);
					}
					else if (strcmp(words[0], "Kai") == 0)
					{
						lastMatParas->Kai = StringToDouble(words[2]);
					}
					else if (strcmp(words[0], "SpecificHeat") == 0)
					{
						lastMatParas->SpecificHeat = StringToDouble(words[2]);
						lastMatParas->SpecificHeat *= 1e-5; // Convert the unit of the input file J/(g*K) to the standard unit in the program 10^5 J/(g*K)
					}
					else if (strcmp(words[0], "Gamma") == 0)
					{
						lastMatParas->Gamma = StringToDouble(words[2]);
					}
					else if (strcmp(words[0], "MolMass") == 0)
					{
						lastMatParas->MolMass = StringToDouble(words[2]);
					}
					else if (strcmp(words[0], "MotionViscos") == 0)
					{
						lastMatParas->MotionViscos = StringToDouble(words[2]);
						lastMatParas->MotionViscos *= 0.01; // Convert the unit m^2/s of the input file to the standard unit cm^2/us in the program
					}
				}
			}
		}
	}
}


/// Converts the material number entered by the user in the input file to the system material position number.
int ConvertUserInputMatNumToSystemMatID(int matNum)
{
	int systemMatID;
	for (systemMatID = matNum; systemMatID < MatParasListLen; systemMatID++)
	{
		if (MatParasList[systemMatID].MaterialNumber == matNum) break;
	}
	return (systemMatID);
}


/// Equation of state, calculation of pressure and specific internal energy from density, temperature
void EOSCalcPressureEnergy(int matId, double density, double temperature, double * pressure, double * specificEnergy)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	EOSTypeEnum eosType = material->EOSType;
	//////////////////////////////////////////// Equation of state for material 1 (positive fluid equation of state)
	if (eosType == EOSTypeEnum_Linear)
	{
		(*pressure) = material->Kai * (density - material->NormalDensity) / material->NormalDensity;
		(*specificEnergy) = temperature * material->SpecificHeat;
	}
	//////////////////////////////////////////// Equation of state for material 2 (polygonal gas equation of state)
	else if (eosType == EOSTypeEnum_IdealGas)
	{
		(*pressure) = temperature * ROfGas * density / material->MolMass;
		(*specificEnergy) = (*pressure) / ((material->Gamma - 1) * density);
	}
	////////////////////////////////////////////
	else
	{
		(*pressure) = 0.0;
		(*specificEnergy) = 0.0;
	}
}

/// Equation of state, calculation of pressure and temperature from density, specific internal energy
void EOSCalcPressureTemperature(int matId, double density, double specificEnergy, double * pressure, double * temperature)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	EOSTypeEnum eosType = material->EOSType;
    //////////////////////////////////////////// Equation of state for material 1 (positive fluid equation of state)
    if (eosType == EOSTypeEnum_Linear)
    {
        (*pressure) = material->Kai * (density - material->NormalDensity) / material->NormalDensity;
        (*temperature) = specificEnergy / material->SpecificHeat;
    }
    //////////////////////////////////////////// Equation of state for material 2 (polygonal gas equation of state)
    else if (eosType == EOSTypeEnum_IdealGas)
    {
        (*pressure) = specificEnergy * ((material->Gamma - 1) * density);
        (*temperature) = (*pressure) * material->MolMass / (ROfGas * density);
    }
    ////////////////////////////////////////////
    else
    {
        (*pressure) = 0.0;
        (*temperature) = 0.0;
    }
}


/// Using the dichotomy method, calculation density from pressure and temperature
double EOSCalcDensity(int matId, double pressure, double temperature)
{
	double pTemp, specificEnergy, density, d1, d2;
	int i;
	///////////////////////////////////////////////////// First get the upper and lower bounds of the density
	density = 1.0;
	EOSCalcPressureEnergy(matId, density, temperature, &pTemp, &specificEnergy);
	if (pTemp < pressure)
	{
		d1 = density;
		//////////////
		for (i = 0; i < 100; i++)
		{
			d2 = d1 * 2;
			EOSCalcPressureEnergy(matId, d2, temperature, &pTemp, &specificEnergy);
			if (pTemp > pressure) break;
			d1 = d2;
		}
	}
	else
	{
		d2 = density;
		//////////////
		for (i = 0; i < 20; i++)
		{
			d1 = d2 / 2;
			EOSCalcPressureEnergy(matId, d1, temperature, &pTemp, &specificEnergy);
			if (pTemp < pressure) break;
			d2 = d1;
		}
	}
	///////////////////////////////////////////////////// Then search with dichotomy until pTemp is very close to pressure
	do
	{
		density = 0.5 * (d1 + d2);
		EOSCalcPressureEnergy(matId, density, temperature, &pTemp, &specificEnergy);
		if (fabs(pTemp - pressure) < 0.001 * pressure || d2 - d1 < 0.0000001 * d2) break;
		////////////////////////////////////
		if (pTemp > pressure) d2 = density;
		else d1 = density;
	} while (true);
	//////////////////////////
	return (density);
}


/// Computes viscous stress from the material's strain rate and density.
double CalcViscosStress(int matId, double strainRatio, double density, double cViscDynArti)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	double viscosStress = 0.0;
	////////////////////////////////////////////
	if (material->MotionViscos > cViscDynArti)
	{
		// Calculating viscous stress from kinematic viscosity coefficient requires multiplying material density
		viscosStress = -material->MotionViscos * strainRatio * density; 
	}
	else
	{
		viscosStress = -cViscDynArti * strainRatio * density;
	}
	////////////////////////////////////////////
	return (viscosStress);
}


/// Returns the viscosity coefficient of the material
double GetViscousCoefficient(int matId)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	return (material->MotionViscos);
}


/// Calculate the speed of sound for a given material number from pressure and specific mass density
double CalcSoundVelocity(int matId, double density, double temperature)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	double soundVelocity = 1.0;
	/////////////////////////////////////////////////
	EOSTypeEnum eosType = material->EOSType;
	//////////////////////////////////////////// 
	if (eosType == EOSTypeEnum_Linear)
	{
		soundVelocity = sqrt(material->Kai / material->NormalDensity);
	}
	//////////////////////////////////////////// 
	else if (eosType == EOSTypeEnum_IdealGas)
	{
		soundVelocity = sqrt(material->Gamma * temperature * ROfGas / material->MolMass);
	}
	////////////////////////////////////////////
	return (soundVelocity);
}


/// Calculate de/dT for a given material
double CalcSpecificEnergyToTemperature(int matId, double density, double temperature)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	double deOverDT = 1.0;
	/////////////////////////////////////////////////
	EOSTypeEnum eosType = material->EOSType;
	//////////////////////////////////////////// 
	if (eosType == EOSTypeEnum_Linear)
	{
		deOverDT = material->SpecificHeat;
	}
	//////////////////////////////////////////// 
	else if (eosType == EOSTypeEnum_IdealGas)
	{
		deOverDT = ROfGas / ((material->Gamma - 1) * material->MolMass);
	}
	////////////////////////////////////////////
	return (deOverDT);
}


/// Calculate the viscous stress tensor from the material's strain rate and density
tensor2D CalcViscosStressTensor(int matId, tensor2D strainRatioTensor, double density, double cViscDynArti)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	tensor2D viscosStress;
	//////////////////////////////////////////// 
	if (matId == 0)
	{
		// Calculating viscous stress from kinematic viscosity coefficient requires multiplying material density
		viscosStress = CValueMultTensor2D(-cViscDynArti * density, strainRatioTensor); 
	}
	else
	{
		// Physical viscosity and artificial dynamic viscosity take the maximum value
		double maxViscCoeff = max(material->MotionViscos, cViscDynArti);
		// scalar viscosity
		viscosStress = CValueToTensor2D((-maxViscCoeff * density) * TraceAverage(strainRatioTensor));	
	}
	////////////////////////////////////////////
	return (viscosStress);
}

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



// double TableEOSCalcDensity(int eosTableIndex, double pressure, double temperature);
/// <summary>
/// 创建真空材料
/// </summary>
/// <returns></returns>
void MakeVacuumMaterial(Material matVacuum)
{
	/////////////////////////////////////////////
	matVacuum->MaterialNumber = 0;
	strcpy(matVacuum->MaterialName, "Vacuum");
	///////////////////////
	matVacuum->EOSType = EOSTypeEnum_Linear;
	matVacuum->NormalDensity = 1e-8;    //室温室压密度
	matVacuum->Kai = 1e-10;             //体模量（可以在计算中动态调整）
	matVacuum->SpecificHeat = 1.0;
}


/// <summary>
/// 从文件中导入材料数据列表
/// </summary>
/// <param name="sr"></param>
void LoadMaterialListFromFile(FILE * fr)
{
	char line[1001], words[3][1001];
	////////////////////////////////////////////// 第 0 个元素放真空材料
	MakeVacuumMaterial(&MatParasList[MatParasListLen++]);
	////////////////////////////////////////////// 找到“>>”记号
	if (SearchNextInFile(fr, ">> 材料设定") == false)
	{
		ERR("没有找到“>> 材料设定”\n");
	}
	////////////////////////////////////////////// 读取材料列表
	{
		Material lastMatParas = NULL;
		int materialNumber = 1;
		bool relateMaterialIsNeed = false;//标记是否需要创建关联材料
		while (true)
		{
			///////////////////////////////////////// 读取一个材料
			fgets(line, 1000, fr);
			if (feof(fr))
			{
				ERR("材料格式不符合要求\n");
			}
			if (ContainStr(line, "--------------------"))
			{
				////////////////////////////////////////// 创建一个新材料
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
						lastMatParas->SpecificHeat *= 1e-5;//将输入文件的单位 J/(g*K) 转换为程序中的标准单位 10^5 J/(g*K)
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
						lastMatParas->MotionViscos *= 0.01;//将输入文件的单位 m^2/s 转换为程序中的标准单位 cm^2/us
					}
				}
			}
		}
	}
}


/// <summary>
/// 将输入文件中用户输入的材料号转换为系统材料位置号。
/// </summary>
int ConvertUserInputMatNumToSystemMatID(int matNum)
{
	int systemMatID;
	for (systemMatID = matNum; systemMatID < MatParasListLen; systemMatID++)
	{
		if (MatParasList[systemMatID].MaterialNumber == matNum) break;
	}
	return (systemMatID);
}


/// <summary>
/// 状态方程，从密度、温度计算压强和比内能
/// </summary>
/// <param name="material">材料号</param>
/// <param name="density">密度(g/cm^3)</param>
/// <param name="temperature">温度(K)</param>
/// <param name="pressure">压强(10^11 Pa, Mbar)</param>
/// <param name="specificEnergy">单位质量的内能（10^5 J/g）</param>
/// <returns></returns>
void EOSCalcPressureEnergy(int matId, double density, double temperature, double * pressure, double * specificEnergy)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	EOSTypeEnum eosType = material->EOSType;
	//////////////////////////////////////////// 第 1 种系统材料的状态方程（正压流体状态方程）
	if (eosType == EOSTypeEnum_Linear)
	{
		(*pressure) = material->Kai * (density - material->NormalDensity) / material->NormalDensity;
		(*specificEnergy) = temperature * material->SpecificHeat;
	}
	//////////////////////////////////////////// 第 2 种系统材料的状态方程（多方气体状态方程）
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

/// <summary>
/// 状态方程，从密度、比内能计算压强和温度
/// </summary>
/// <param name="material">材料号</param>
/// <param name="density">密度(g/cm^3)</param>
/// <param name="temperature">温度(K)</param>
/// <param name="pressure">压强(10^11 Pa, Mbar)</param>
/// <param name="specificEnergy">单位质量的内能（10^5 J/g）</param>
/// <returns></returns>
void EOSCalcPressureTemperature(int matId, double density, double specificEnergy, double * pressure, double * temperature)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	EOSTypeEnum eosType = material->EOSType;
    //////////////////////////////////////////// 第 1 种系统材料的状态方程（正压流体状态方程）
    if (eosType == EOSTypeEnum_Linear)
    {
        (*pressure) = material->Kai * (density - material->NormalDensity) / material->NormalDensity;
        (*temperature) = specificEnergy / material->SpecificHeat;
    }
    //////////////////////////////////////////// 第 2 种系统材料的状态方程（多方气体状态方程）
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


/// <summary>
/// 利用二分法，从压强和温度，倒推密度。
/// </summary>
/// <param name="material"></param>
/// <param name="pressure"></param>
/// <param name="temperature"></param>
/// <param name="density"></param>
/// <param name="specificEnergy"></param>
double EOSCalcDensity(int matId, double pressure, double temperature)
{
	double pTemp, specificEnergy, density, d1, d2;
	int i;
	///////////////////////////////////////////////////// 首先获得密度的上下限
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
	///////////////////////////////////////////////////// 然后用二分法搜寻，直到pTemp与pressure相隔很近
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


/// <summary>
/// 从材料的应变率和密度计算粘性应力.
/// 人工动态粘性系数 和 物理粘性系数 中，哪一个更大就用哪一个。
/// </summary>
/// <param name="idMaterial"></param>
/// <param name="strainRatio"></param>
/// <param name="density"></param>
/// <returns></returns>
double CalcViscosStress(int matId, double strainRatio, double density, double cViscDynArti)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	double viscosStress = 0.0;
	////////////////////////////////////////////
	if (material->MotionViscos > cViscDynArti)
	{
		viscosStress = -material->MotionViscos * strainRatio * density;//从运动粘性系数计算粘性应力需要乘以材料密度
	}
	else
	{
		viscosStress = -cViscDynArti * strainRatio * density;
	}
	////////////////////////////////////////////
	return (viscosStress);
}


/// <summary>
/// 读取一种材料的粘性系数。（用于计算粘性决定的时间步长）
/// 20171104
/// </summary>
/// <param name="idMaterial"></param>
/// <returns></returns>
double GetViscousCoefficient(int matId)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	return (material->MotionViscos);
}


/// <summary>
/// 从压强和比质量密度计算指定材料号的声速
/// </summary>
/// <param name="idMaterial"></param>
/// <param name="temperature"></param>
/// <param name="density"></param>
/// <returns></returns>
double CalcSoundVelocity(int matId, double density, double temperature)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	double soundVelocity = 1.0;
	/////////////////////////////////////////////////
	EOSTypeEnum eosType = material->EOSType;
	//////////////////////////////////////////// 第 1 种系统材料的声速（轻重流体状态方程）
	if (eosType == EOSTypeEnum_Linear)
	{
		soundVelocity = sqrt(material->Kai / material->NormalDensity);
	}
	//////////////////////////////////////////// 第 2 种系统材料的声速（多方气体状态方程）。改正声速公式错误20190325（错误源于一维程序）
	else if (eosType == EOSTypeEnum_IdealGas)
	{
		soundVelocity = sqrt(material->Gamma * temperature * ROfGas / material->MolMass);
	}
	////////////////////////////////////////////
	return (soundVelocity);
}


/// <summary>
/// 计算指定材料的de/dT
/// </summary>
/// <param name="idMaterial"></param>
/// <param name="density"></param>
/// <param name="temperature"></param>
/// <returns></returns>
double CalcSpecificEnergyToTemperature(int matId, double density, double temperature)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	double deOverDT = 1.0;
	/////////////////////////////////////////////////
	EOSTypeEnum eosType = material->EOSType;
	//////////////////////////////////////////// 第 1 种系统材料（线弹性流体状态方程）
	if (eosType == EOSTypeEnum_Linear)
	{
		deOverDT = material->SpecificHeat;
	}
	//////////////////////////////////////////// 第 2 种系统材料（多方气体状态方程）
	else if (eosType == EOSTypeEnum_IdealGas)
	{
		deOverDT = ROfGas / ((material->Gamma - 1) * material->MolMass);
	}
	////////////////////////////////////////////
	return (deOverDT);
}


/// <summary>
/// 从材料的应变率和密度计算粘性应力张量
/// </summary>
/// <param name="idMaterial"></param>
/// <param name="strainRatio"></param>
/// <param name="density"></param>
/// <returns></returns>
tensor2D CalcViscosStressTensor(int matId, tensor2D strainRatioTensor, double density, double cViscDynArti)
{
	Material material = matId >= 0 ? &MatParasList[matId] : &MatWall[-matId];
	/////////////////////////////////////////////
	tensor2D viscosStress;
	//////////////////////////////////////////// 从运动粘性系数计算粘性应力需要乘以材料密度
	if (matId == 0)
	{
		viscosStress = CValueMultTensor2D(-cViscDynArti * density, strainRatioTensor);//从运动粘性系数计算粘性应力需要乘以材料密度
	}
	else
	{
		double maxViscCoeff = max(material->MotionViscos, cViscDynArti);//物理粘性和人工动态粘性哪一个更大就用哪一个20180403	
		// 标量粘性
		viscosStress = CValueToTensor2D((-maxViscCoeff * density) * TraceAverage(strainRatioTensor));	
	}
	////////////////////////////////////////////
	return (viscosStress);
}

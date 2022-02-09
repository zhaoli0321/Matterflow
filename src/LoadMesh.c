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
/// <summary>
/// 从输入文件 nodeListFile 读取节点数据。
/// 初始作者：谢龙20150527
/// </summary>
/// <param name="nodeListFile">节点数据所在的文件名</param>
/// <returns></returns>
void LoadNodes(char nodeListFile[], struct nodeFromFile nodesArray[], int * nodesArrLenOut)
{
	char line[1001], words[5][101];
    int nodesArrLen = 0;
    ////////////////////////////////////////////////////////// 1)打开文件
    FILE * fr = MyOpenFile(nodeListFile, "r");
	////////////////////////////////////////////////////////// 2)跳过首行
	if(SearchNextInFile(fr, "速度y") == false)
	{
		ERR("读取格点文件首行出错\n");
	}
    ////////////////////////////////////////////////////////// 2)读取数据
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
    ////////////////////////////////////////////////////////// 3)保存节点数据
    if (nodesArrLen == 0)
    {
        ERR("错误: 文件中不包含有效的节点数据! \n");
    }
    //////////////////////////////////////////////////////////
    fclose(fr);
	(*nodesArrLenOut) = nodesArrLen;
}
/// <summary>
/// 从输入文件 cellListFile 读取三角形单元数据
/// 初始作者：谢龙20150527
/// 增加输入多样性：有无磁场 Bz(Tesla)都兼容，赵梨
/// </summary>
/// <param name="cellListFile">三角形单元数据所在的文件名</param>
/// <param name="nodes">预先准备好的格点集合</param>
/// <param name="phaseVecs">预先准备好的材料（态矢量）集合</param>
/// <returns></returns>
void LoadCells(char cellListFile[], struct triangleFromFile trgsArray[], int * trgsArrLenOut)
{
	char line[1001], words[8][101];
    int trgsArrLen = 0;
    bool isHaveBz = false;
    ////////////////////////////////////////////////////////// 1)打开文件
    FILE * fr = MyOpenFile(cellListFile, "r");
	////////////////////////////////////////////////////////// 2)跳过首行
    // if(SearchNextInFile(fr, "Bz(Tesla)") == true)
	// {
	// 	isHaveBz = true;
    //     printf("Bz(Tesla)\n");
	// }else{
    //     ERR("读取单元文件首行出错\n");
    // }

    fgets(line, 1000, fr);
    if (ContainStr(line, "Bz(Tesla)"))
	{
		isHaveBz = true;
        // printf("Bz(Tesla)\n");
	}else if(ContainStr(line, "温度"))
	{
		isHaveBz = false;
        // printf("温度\n");
	}else{
        ERR("读取单元文件首行出错\n");
    }

    // printf("isHaveBz: %d\n", isHaveBz);
    ////////////////////////////////////////////////////////// 2)读取数据
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
                        trg->Pressure    = StringToDouble(words[5]) * 1e-6;//将输入文件的单位 10^5Pa 转换为程序中的标准单位 100GPa
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
                        trg->Pressure    = StringToDouble(words[5]) * 1e-6;//将输入文件的单位 10^5Pa 转换为程序中的标准单位 100GPa
                        trg->Temperature = StringToDouble(words[6]);
                    }
                }
            }
        }
    }
    ////////////////////////////////////////////////////////// 3)保存单元数据
    if (trgsArrLen == 0)
    {
        ERR("错误: 文件中不包含有效的单元数据!\n");
    }
    //////////////////////////////////////////////////////////
    fclose(fr);
    (*trgsArrLenOut) = trgsArrLen;
}
/// <summary>
/// 对边界处的三角形连接关系作“周期性”处理。
/// 做法上，就是修改边界三角形的顶点。
/// 该函数需要Ansys生成的网格文件满足一定的要求：
///   1. 必须是四方形的区域
///   2. 四方形对边的格点必须数目相等，位置匹配。
/// </summary>
/// <param name="nodes"></param>
/// <param name="trgs"></param>
void  MakeCycleTrgsOnEdge(struct nodeFromFile nodes[], int totNodes, struct triangleFromFile trgs[], int totTrgs)
{
    double minX, maxX, minY, maxY, averageSideLength, tinyLength;
    int * idForEdgeNodes;
    int totEdgeNodes;
	int i,t,k;

    //////////////////////////////////// 先获得网格的最小和最大X，Y坐标以及单元平均边长，便于下面判断格点在哪条边界上。
    minX = maxX = nodes[0].Pos.X;
    minY = maxY = nodes[0].Pos.Y;
    for (i = 0; i < totNodes; i++ )
    {       
        NodeFromFile node = &nodes[i];
        ////////////////////////
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
    averageSideLength = sqrt((maxX - minX) * (maxY - minY) / totTrgs);
    tinyLength = 0.1 * averageSideLength;
    //////////////////////////////////// 然后找出所有边界上的格点备用（主要是为了加快下文的格点遍历速度）
    idForEdgeNodes = (int *)malloc(sizeof(int)*totNodes);//new int[totNodes];
	totEdgeNodes = 0;
    for (i = 0; i < totNodes; i++)
    {
        NodeFromFile node = &nodes[i];
        ////////////////////////
        double posX = node->Pos.X;
        double posY = node->Pos.Y;
        if (fabs(posX - minX) < tinyLength || fabs(posX - maxX) < tinyLength
            || fabs(posY - minY) < tinyLength || fabs(posY - maxY) < tinyLength)
        {
            idForEdgeNodes[totEdgeNodes++] = i;
        }
    }
    //////////////////////////////////// 然后修改边界三角形的顶点属性（我们这里选择把右边和上边边界上的顶点废弃掉）
    for (t = 0; t < totTrgs; t++ )
    {
        TriangleFromFile trg = &trgs[t];
        ////////////////////////////////
        for (i = 0; i < 3; i++)
        {
            int newID = trg->NodeIDs[i];           
            //////////////////////////// 该顶点如果在右边边界上，就改为指向左边界的对应格点。
            if (fabs(nodes[newID].Pos.X - maxX) < tinyLength)
            {
                for (k = 0; k < totEdgeNodes; k++ )
                {
                    int idEdge = idForEdgeNodes[k];
                    if (fabs(nodes[idEdge].Pos.Y - nodes[newID].Pos.Y) < tinyLength && fabs(nodes[idEdge].Pos.X - minX) < tinyLength)
                    {
                        newID = idEdge;
                        break;
                    }
                }
            }
            //////////////////////////// 该顶点如果在上边边界上，就改为指向下边界的对应格点。
            if (fabs(nodes[newID].Pos.Y - maxY) < tinyLength)
            {
                for (k = 0; k < totEdgeNodes; k++)
                {
                    int idEdge = idForEdgeNodes[k];
                    if (fabs(nodes[idEdge].Pos.X - nodes[newID].Pos.X) < tinyLength && fabs(nodes[idEdge].Pos.Y - minY) < tinyLength)
                    {
                        newID = idEdge;
                        break;
                    }
                }
            }
            ////////////////////////////
            trg->NodeIDs[i] = newID;
        }
    }
	//////////////////////////////////////////////////////////////// free
	free(idForEdgeNodes);
}

/// <summary>
/// 从格点和三角单元集合构建三角形邻域关系。
/// （算法是钟敏和谢龙教给我的。20151109）
/// </summary>
/// <param name="nodes"></param>
/// <param name="trgs"></param>
void ConstructTrgNeighbours(struct nodeFromFile nodes[], int nodesLen, struct triangleFromFile trgs[], int trgsLen)
{
	int i, j, k, t;
    //////////////////////////////////////////////////////////// 第一步,构建每个节点的邻域三角形
	/// 先把格点的邻域三角形数目清零
	for (i = 0; i < nodesLen; i++)
	{
		nodes[i].TrgIDsLen = 0;
	}
    /// 遍历三角单元，构建每个节点的邻域三角形
    for (i = 0; i < trgsLen; i++)
    {
        int * nodeIds = trgs[i].NodeIDs;
        nodes[nodeIds[0]].TrgIDs[nodes[nodeIds[0]].TrgIDsLen++] = i;
        nodes[nodeIds[1]].TrgIDs[nodes[nodeIds[1]].TrgIDsLen++] = i;
        nodes[nodeIds[2]].TrgIDs[nodes[nodeIds[2]].TrgIDsLen++] = i;
    }
    //////////////////////////////////////////////////////////// 第二步
    /// 再次遍历所有三角单元，从每条边的两个节点的共有三角形判断
    /// 该三角形这条边的邻域三角形。
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
/// <summary>
/// 从Ansys格点和单元生成程序中所需的标准格点和单元.
/// 20160424修正：
///   发现Ansys生成的格点和三角单元文件中，格点文件中的格点有许多多余的。
///   在这个函数中，将消除这些多余格点。
/// </summary>
/// <param name="nodesFromFile"></param>
/// <param name="trgsFromFile"></param>
/// <param name="verts"></param>
/// <param name="trgs"></param>
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
    /// 三角单元格点的设置
    for (i = 0; i < totTrgs; i++)
    {
        int* nodeIDs = trgsFromFile[i].NodeIDs;
        ////////////////////////////////////////
        trgs[i].Vertices[0] = &verts[nodeIDs[0]];
        trgs[i].Vertices[1] = &verts[nodeIDs[1]];
        trgs[i].Vertices[2] = &verts[nodeIDs[2]];
        ////////////////////////////////////////       
    }

    ///////////////////////////////////////////////////////////// 设置三角形的物理属性
    /// 三角形质量、内能和标准边长
    int* neighbourIDs;
    for (i = 0; i < totTrgs; i++)
    {
        Triangle trg = &trgs[i];

        /// 创建邻边三角形单元
        // printf("trg[%d] = ", i);
        neighbourIDs = trgsFromFile[i].NeighbourTrgIDs;
        for (j = 0; j < 3; j++)
        {
            if (neighbourIDs[j] != -1)
            {
                trgs[i].NeighbourTrgs[j] = &trgs[neighbourIDs[j]];
                // printf("%d ", trgs[i].NeighbourTrgs[j]->Index);
            }
            else
            {
                trgs[i].NeighbourTrgs[j] = NullTriangle;
                // if (trgs[i].NeighbourTrgs[j] == NullTriangle) printf("NULL ");
            }
        }
        // printf("\n");

        ////////////////////////////////////////////// 计算三角单元面积
		double area;
        area = CalcTrgArea(trg->Vertices[0]->Pos, trg->Vertices[1]->Pos, trg->Vertices[2]->Pos);
		
        ////////////////////////////////////////////// 设置密度、温度、磁通量这些基本物理量和其它物理量
        if (trgsFromFile[i].MaterialId > 0)
        {
            ////////////////////////////// 基本物理量（密度、温度、电流密度）
            Material material = &MatParasList[trgsFromFile[i].MaterialId];
            {
                trg->MaterialId = trgsFromFile[i].MaterialId;
                //////////////////////////////////////////////
                trg->Temperature = trgsFromFile[i].Temperature;
                trg->Density = EOSCalcDensity(trg->MaterialId, trgsFromFile[i].Pressure, trgsFromFile[i].Temperature);
            }
            ////////////////////////////// 其它物理量
            trg->Mass = trg->Density * area;
            {
                double pressure, specificEnergy;
                EOSCalcPressureEnergy(trg->MaterialId, trg->Density, trg->Temperature, &pressure, &specificEnergy);
                trg->InternalEnergy = trg->Mass * specificEnergy;
            }
        }
        else if (trgsFromFile[i].MaterialId == 0)//如果是真空单元，就采用密度极小的线弹性流体描述
        {
            trg->MaterialId = 0;
            trg->Density = MatParasList[0].NormalDensity;
            trg->Mass = trg->Density * area;
        }
        else//如果是“墙壁材料”单元，也采用线弹性流体描述。（其实从逻辑上来说，墙壁单元不需要质量和密度属性，但是为了网格重分重映代码的简洁一致性，还是给一个大于零的值）
        {
            trg->MaterialId = trgsFromFile[i].MaterialId;
            trg->Density = MatWall[-trg->MaterialId].NormalDensity;
            trg->Mass = trg->Density * area;
		}
    }
    ///////////////////////////////////////////////////////////// 然后设置格点的物理属性
    /// 格点质量
    for (i = 0; i < totTrgs; i++)
    {
        Triangle trg = &trgs[i];
        if (trg->IsDead) continue;
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

/// <summary>
/// 建立边界格点的 OnUpDownBoundary 和 OnLeftRightBoundary 属性
/// </summary>
/// <param name="vertices"></param>
void SetUpDownLeftRightBoundaryVertices(struct vertex vertices[], int vertsArrLen)
{
	int i;
    //////////////////////////////////// 先获得 minX, minY, tinyLength ，便于下面判断格点是否在边界上。
    /// 注意，这里配合了上文：（我们这里选择把右边和上边边界上的顶点废弃掉）
    /// 因此只需要用到 minX, minY 来判断是否在边界上
    double minX, minY, maxX, maxY, tinyLength;
    {
        minX = maxX = vertices[0].Pos.X;
        minY = maxY = vertices[0].Pos.Y;
        for (i = 0; i < vertsArrLen; i++ )
        {
            Vertex vert = &vertices[i];
            if (vert->IsDead) continue;
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
    //////////////////////////////////// 下面判断格点是否在边界上。
    for (i = 0; i < vertsArrLen; i++)
    {
        Vertex vert = &vertices[i];
        if (vert->IsDead) continue;
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
/// <summary>
/// 将坐标平移——左下角归零。并借此机会确定网格物体宽度和高度.
/// 然后，将网格沿轴线（r=0）镜像翻倍。（注意，轴心处的格点也被翻倍了，但没有三角形指向这些多出的格点，所以将在后面的操作中被去除掉。）
/// </summary>
/// <param name="nodes"></param>
/// <param name="width"></param>
/// <param name="height"></param>
void TranslatePositions(struct nodeFromFile nodes[], int totNodes, struct triangleFromFile trgs[], int totTrgs, double * widthOut, double * heightOut)
{
	int i;
    ///////////////////////////////////////////////////////////////////// 第一步：将坐标平移——左下角归零。并借此机会确定网格物体宽度和高度.
    /// 根据最大、最小坐标值确定平移量
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
    //////////////////////// 坐标平移(左下角归零)
	{
		vec2D translation = MakeVec2D(-minX, -minY);
		for ( i = 0; i < totNodes; i++)
		{
			NodeFromFile node = &nodes[i];
			//坐标的平移
			Vec2DSelfAdd(&node->Pos, translation);
		}
	}
    //////////////////////// 区域宽度（半宽）和高度
    (*widthOut) = maxX - minX;
    (*heightOut) = maxY - minY;
}
/// <summary>
/// 从Ansys生成的节点和网格文件中导入物体
/// </summary>
/// <param name="nodesFileName"></param>
/// <param name="cellsFileName"></param>
/// <returns></returns>
void SetMeshFromMeshFiles(char nodesFileName[], char cellsFileName[])
{
    ///////////////////////////////////////////////// A.1 导入Ansys节点和单元数据
    struct nodeFromFile * nodesFromFile = (struct nodeFromFile *)malloc(sizeof(struct nodeFromFile) * MAX_VERTS);//new NodeFromFile[MeshObject.MAX_VERTS];
    struct triangleFromFile * trgsFromFile = (struct triangleFromFile *)malloc(sizeof(struct triangleFromFile) * MAX_VERTS * 2);////new TriangleFromFile[MeshObject.MAX_VERTS * 2];
    int totNodes, totTrgs;
    LoadNodes(nodesFileName, nodesFromFile, &totNodes);
    LoadCells(cellsFileName, trgsFromFile, &totTrgs);
    ///////////////////////////////////////////////// A.2 将坐标平移——左下角归零。（并借此机会确定网格物体宽度和高度）。然后网格镜像翻倍。
	TranslatePositions(nodesFromFile, totNodes, trgsFromFile, totTrgs, &MeshObj.Width, &MeshObj.Height);
    ///////////////////////////////////////////////// A.3 构建三角形邻域关系
    ConstructTrgNeighbours(nodesFromFile, totNodes, trgsFromFile, totTrgs);
    ///////////////////////////////////////////////// B.1 从建模文件格点和单元生成程序中所需的标准格点和单元
    ConvertModelMeshToStandard(nodesFromFile, totNodes, trgsFromFile, totTrgs);
    ///////////////////////////////////////////////// B.2 建立边界格点的 OnUpDownBoundary 和 OnLeftRightBoundary 属性
    SetUpDownLeftRightBoundaryVertices(MeshObj.Vertices, MeshObj.VertsArrLen);
    ///////////////////////////////////////////////// 设置所有因变量
    // SetAllDependentVariablesOfTrgs();
    /////////////////////////////////////////////////
	free(nodesFromFile);
	free(trgsFromFile);
}
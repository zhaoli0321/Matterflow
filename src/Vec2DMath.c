/*! \file  Vec2DMath.c
 *
 *  \brief Define some functions for Vec2D
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"

vec2D MakeVec2D(double r, double z)
{
	vec2D out;
	out.X = r;
	out.Y = z;
	return(out);
}

void SelfRotate90(vec2D * v)
{
	double temp = v->X;
	v->X = -v->Y;
	v->Y = temp;
}

vec2D Rotate90(vec2D v)
{
	vec2D vec;
	vec.X = -v.Y;
	vec.Y = v.X;
	return (vec);
}

vec2D Rotate270(vec2D v)
{
	vec2D vec;
	vec.X = v.Y;
	vec.Y = -v.X;
	return (vec);
}

/// <summary>
/// 
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
vec2D Vec2DAdd(vec2D v1, vec2D v2)
{
	v1.X += v2.X;
	v1.Y += v2.Y;
	return (v1);
}

/// <summary>
/// 
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
vec2D Vec2DAdd3(vec2D v1, vec2D v2, vec2D v3)
{
	v1.X += v2.X + v3.X;
	v1.Y += v2.Y + v3.Y;
	return (v1);
}

/// <summary>
/// 
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
void Vec2DSelfAdd(vec2D * v1, vec2D v2)
{
	(*v1).X += v2.X;
	(*v1).Y += v2.Y;
}

/// <summary>
/// 
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
vec2D Vec2DSub(vec2D v1, vec2D v2)
{
	v1.X -= v2.X;
	v1.Y -= v2.Y;
	return (v1);  // 因为结构体是按值传参数的，所以可以直接借用形参v1。
}

/// <summary>
/// 
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
void Vec2DSelfSub(vec2D * v1, vec2D v2)
{
	(*v1).X -= v2.X;
	(*v1).Y -= v2.Y;
}

/// <summary>
/// 两个位置点之间的距离
/// </summary>
/// <param name="pos1"></param>
/// <param name="pos2"></param>
/// <returns></returns>
double Dist(vec2D pos1, vec2D pos2)
{
	return (sqrt((pos1.X - pos2.X) * (pos1.X - pos2.X) + (pos1.Y - pos2.Y) * (pos1.Y - pos2.Y)));
}

/// <summary>
/// 
/// </summary>
/// <param name="vec"></param>
/// <returns></returns>
double CalcLength(vec2D vec)
{
	return (sqrt(vec.X * vec.X + vec.Y * vec.Y));
}

double CalcLengthSquar(vec2D vec)
{
	return (vec.X * vec.X + vec.Y * vec.Y);
}

/// <summary>
/// 
/// </summary>
/// <param name="vec"></param>
/// <param name="ratio"></param>
/// <returns></returns>
vec2D Vec2DMultCValue(vec2D vec, double ratio)
{
	vec.X *= ratio;
	vec.Y *= ratio;
	return (vec);  // 因为结构体是按值传参数的，所以可以直接借用形参vec。
}

/// <summary>
/// 
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
void Vec2DSelfMultCValue(vec2D * v1, double c)
{
	(*v1).X *= c;
	(*v1).Y *= c;
}

/// <summary>
/// 
/// </summary>
/// <param name="ratio"></param>
/// <param name="vec"></param>
/// <returns></returns>
vec2D CValueMultVec2D(double ratio, vec2D vec)
{
	vec.X *= ratio;
	vec.Y *= ratio;
	return (vec);  // 因为结构体是按值传参数的，所以可以直接借用形参vec。
}

/// <summary>
/// 二维矢量叉乘。（注意，二维矢量叉乘返回的是一个实数，而不是矢量。）
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
double Cross(vec2D v1, vec2D v2)
{
	return (v1.X * v2.Y - v1.Y * v2.X);
}

/// <summary>
/// 二维矢量点乘。
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
double Dot(vec2D v1, vec2D v2)
{
	return (v1.X * v2.X + v1.Y * v2.Y);
}

/// <summary>
/// 计算三角形的面积
/// 这是一个有正负的面积，当顶点1,2,3为逆时针排布时，面积为正。
/// 可以用面积的正负来判断网格有没有发生畸变。
/// </summary>
/// <param name="pos0"></param>
/// <param name="pos1"></param>
/// <param name="pos2"></param>
/// <returns></returns>
double CalcTrgArea(vec2D pos0, vec2D pos1, vec2D pos2)
{
	return (0.5 * Cross(Vec2DSub(pos1, pos0), Vec2DSub(pos2, pos0)));
}

/// <summary>
/// 
/// </summary>
/// <param name="vec"></param>
/// <param name="ratio"></param>
/// <returns></returns>
vec2D Vec2DDivideCValue(vec2D vec, double ratio)
{
	vec.X /= ratio;
	vec.Y /= ratio;
	return (vec);  // 因为结构体是按值传参数的，所以可以直接借用形参vec。
}

/// <summary>
/// 从v1指向v2的单位方向矢量
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
vec2D DirectVec(vec2D v1, vec2D v2)
{
	vec2D vec;
	vec = Vec2DSub(v2, v1);
	{
		double ratio = CalcLength(vec);
		vec.X /= ratio;
		vec.Y /= ratio;
	}
	return (vec);
}

/// <summary>
/// 矢量的单位长度方向矢量
/// </summary>
/// <param name="vec"></param>
/// <returns></returns>
vec2D UnitVec(vec2D vec)
{
	double ratio = CalcLength(vec);
	vec.X /= ratio;
	vec.Y /= ratio;
	return (vec);  // 因为结构体是按值传参数的，所以可以直接借用形参vec。
}

/// <summary>
/// 计算三角形的三条高中的高度相对变化率最大的值。
/// </summary>
/// <param name="r1"></param>
/// <param name="r2"></param>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
double CalcTrgMaxHeightChangeRate(vec2D r1, vec2D r2, vec2D v1, vec2D v2)
{
	////////////////////////////////////// 第一步，算出应变率张量
	tensor2D strainRatio = CalcStrainRateTensor(r1, r2, v1, v2);
	////////////////////////////////////// 第二步，应变率张量对角化后，取最大对角元
	tensor2D diagTensor = Diagonalize(strainRatio);
	double maxRatio = max(fabs(diagTensor.A[0][0]), fabs(diagTensor.A[1][1]));
	//////////////////////////////////////
	return (maxRatio);
}

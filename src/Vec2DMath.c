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

/// Create vector
vec2D MakeVec2D(double r, double z)
{
	vec2D out;
	out.X = r;
	out.Y = z;
	return(out);
}

/// Rotate the vector 90 degree (self)
void SelfRotate90(vec2D * v)
{
	double temp = v->X;
	v->X = -v->Y;
	v->Y = temp;
}

/// Rotate the vector 90 degree
vec2D Rotate90(vec2D v)
{
	vec2D vec;
	vec.X = -v.Y;
	vec.Y = v.X;
	return (vec);
}

/// Rotate the vector 270 degree
vec2D Rotate270(vec2D v)
{
	vec2D vec;
	vec.X = v.Y;
	vec.Y = -v.X;
	return (vec);
}

/// Vector addition
vec2D Vec2DAdd(vec2D v1, vec2D v2)
{
	v1.X += v2.X;
	v1.Y += v2.Y;
	return (v1);
}

/// Three vectors addition
vec2D Vec2DAdd3(vec2D v1, vec2D v2, vec2D v3)
{
	v1.X += v2.X + v3.X;
	v1.Y += v2.Y + v3.Y;
	return (v1);
}

/// Vector self addition
void Vec2DSelfAdd(vec2D * v1, vec2D v2)
{
	(*v1).X += v2.X;
	(*v1).Y += v2.Y;
}

/// Vector subtraction
vec2D Vec2DSub(vec2D v1, vec2D v2)
{
	v1.X -= v2.X;
	v1.Y -= v2.Y;
	return (v1);  
}

/// Vector self subtraction
void Vec2DSelfSub(vec2D * v1, vec2D v2)
{
	(*v1).X -= v2.X;
	(*v1).Y -= v2.Y;
}

/// Distance between two vectors
double Dist(vec2D pos1, vec2D pos2)
{
	return (sqrt((pos1.X - pos2.X) * (pos1.X - pos2.X) + (pos1.Y - pos2.Y) * (pos1.Y - pos2.Y)));
}

/// Modulo length of a vector
double CalcLength(vec2D vec)
{
	return (sqrt(vec.X * vec.X + vec.Y * vec.Y));
}

/// The square of the modulo length of the vector
double CalcLengthSquar(vec2D vec)
{
	return (vec.X * vec.X + vec.Y * vec.Y);
}

/// Vector multiplied by scalar
vec2D Vec2DMultCValue(vec2D vec, double ratio)
{
	vec.X *= ratio;
	vec.Y *= ratio;
	return (vec);  
}

/// Vector multiplied by scalar (self)
void Vec2DSelfMultCValue(vec2D * v1, double c)
{
	(*v1).X *= c;
	(*v1).Y *= c;
}

/// Vector multiplied by scalar
vec2D CValueMultVec2D(double ratio, vec2D vec)
{
	vec.X *= ratio;
	vec.Y *= ratio;
	return (vec);  
}

/// 2D vector cross product. 
/// Note that the 2D vector cross product returns a real number, not a vector.
double Cross(vec2D v1, vec2D v2)
{
	return (v1.X * v2.Y - v1.Y * v2.X);
}

/// 2D vector dot product.
double Dot(vec2D v1, vec2D v2)
{
	return (v1.X * v2.X + v1.Y * v2.Y);
}

/// Calculate the area of a triangle
/// This is a positive and negative area. When vertices 1, 2, and 3 are arranged counterclockwise, the area is positive.
double CalcTrgArea(vec2D pos0, vec2D pos1, vec2D pos2)
{
	return (0.5 * Cross(Vec2DSub(pos1, pos0), Vec2DSub(pos2, pos0)));
}

/// Vector divided by constant
vec2D Vec2DDivideCValue(vec2D vec, double ratio)
{
	vec.X /= ratio;
	vec.Y /= ratio;
	return (vec);  
}

/// Unit direction vector from v1 to v2
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

/// Vector unit length direction vector
vec2D UnitVec(vec2D vec)
{
	double ratio = CalcLength(vec);
	vec.X /= ratio;
	vec.Y /= ratio;
	return (vec);  
}

/// Calculate the maximum value of the relative change rate for the three heights of the triangle.
double CalcTrgMaxHeightChangeRate(vec2D r1, vec2D r2, vec2D v1, vec2D v2)
{
	// The first step is to calculate the strain rate tensor
	tensor2D strainRatio = CalcStrainRateTensor(r1, r2, v1, v2);
	// In the second step, after the strain rate tensor is diagonalized, take the largest diagonal element
	tensor2D diagTensor = Diagonalize(strainRatio);
	double maxRatio = max(fabs(diagTensor.A[0][0]), fabs(diagTensor.A[1][1]));
	//////////////////////////////////////
	return (maxRatio);
}

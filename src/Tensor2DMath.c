/*! \file  Tensor2DMath.c
 *
 *  \brief Define some functions for Tensor2DMath
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"

/// Tensor addition
tensor2D Tensor2DAdd(tensor2D t1, tensor2D t2)
{
    tensor2D t3;
    double (*a1)[2] = t1.A;
    double (*a2)[2] = t2.A;
    double (*a3)[2] = t3.A;
    a3[0][0] = a1[0][0] + a2[0][0];
    a3[0][1] = a1[0][1] + a2[0][1];
    a3[1][0] = a1[1][0] + a2[1][0];
    a3[1][1] = a1[1][1] + a2[1][1];
    return (t3);
}

/// Tensor self addition
void Tensor2DSelfAdd(tensor2D * t1, tensor2D t2)
{
    double (*a1)[2] = t1->A;
    double (*a2)[2] = t2.A;
    a1[0][0] += a2[0][0];
    a1[0][1] += a2[0][1];
    a1[1][0] += a2[1][0];
    a1[1][1] += a2[1][1];
}

/// Tensor subtraction
tensor2D Tensor2DSub(tensor2D t1, tensor2D t2)
{
    tensor2D t3;
    double(*a1)[2] = t1.A;
    double(*a2)[2] = t2.A;
    double(*a3)[2] = t3.A;
    a3[0][0] = a1[0][0] - a2[0][0];
    a3[0][1] = a1[0][1] - a2[0][1];
    a3[1][0] = a1[1][0] - a2[1][0];
    a3[1][1] = a1[1][1] - a2[1][1];
    return (t3);
}

/// Add a real number to a tensor 
/// (essentially multiply a real number by a unit tensor and add it to the tensor)
tensor2D CValueAddTensor2D(double c, tensor2D t2)
{
    tensor2D t3;
    double(*a2)[2] = t2.A;
    double(*a3)[2] = t3.A;
    a3[0][0] = c + a2[0][0];
    a3[0][1] = a2[0][1];
    a3[1][0] = a2[1][0];
    a3[1][1] = c + a2[1][1];
    return (t3);
}

/// Add a real number to a tensor 
/// (essentially multiply a real number by a unit tensor and add it to the tensor)
tensor2D Tensor2DAddCValue(tensor2D t2, double c)
{
    tensor2D t3;
    double(*a2)[2] = t2.A;
    double(*a3)[2] = t3.A;
    a3[0][0] = c + a2[0][0];
    a3[0][1] = a2[0][1];
    a3[1][0] = a2[1][0];
    a3[1][1] = c + a2[1][1];
    return (t3);
}

/// Self add a real number to a tensor 
/// (essentially multiply a real number by a unit tensor and add it to the tensor)
void Tensor2DSelfAddCValue(tensor2D * t1, double c)
{
	double(*a1)[2] = (*t1).A;
	a1[0][0] += c;
	a1[1][1] += c;
}

/// Real multiplied by tensor
tensor2D CValueMultTensor2D(double factor, tensor2D tensor)
{
    tensor2D result;
    double (*a1)[2] = tensor.A, (*a2)[2] = result.A;
    a2[0][0] = a1[0][0] * factor;
    a2[0][1] = a1[0][1] * factor;
    a2[1][0] = a1[1][0] * factor;
    a2[1][1] = a1[1][1] * factor;
    return (result);
}

/// Real multiplied by tensor
tensor2D Tensor2DMultCValue(tensor2D tensor, double factor)
{
    tensor2D result;
    double (*a1)[2] = tensor.A, (*a2)[2] = result.A;
    a2[0][0] = a1[0][0] * factor;
    a2[0][1] = a1[0][1] * factor;
    a2[1][0] = a1[1][0] * factor;
    a2[1][1] = a1[1][1] * factor;
    return (result);
}

/// Real multiplied by tensor (self)
void Tensor2DSelfMultCValue(tensor2D * tensor, double factor)
{
    double (*a)[2] = tensor->A;
    a[0][0] *= factor;
    a[0][1] *= factor;
    a[1][0] *= factor;
    a[1][1] *= factor;
}

/// Tensor multiplied by vector
vec2D Tensor2DMultVec(tensor2D tensor, vec2D vec)
{
    vec2D result;
    double (*a)[2] = tensor.A;
    result.X = a[0][0] * vec.X + a[0][1] * vec.Y;
    result.Y = a[1][0] * vec.X + a[1][1] * vec.Y;
    return (result);
}

/// Tensor multiplied by tensor
tensor2D Tensor2DMult(tensor2D t1, tensor2D t2)
{
    tensor2D t3;
    double(*a1)[2] = t1.A;
    double(*a2)[2] = t2.A;
    double(*a3)[2] = t3.A;
    a3[0][0] = a1[0][0] * a2[0][0] + a1[0][1] * a2[1][0];
    a3[0][1] = a1[0][0] * a2[0][1] + a1[0][1] * a2[1][1];
    a3[1][0] = a1[1][0] * a2[0][0] + a1[1][1] * a2[1][0];
    a3[1][1] = a1[1][0] * a2[0][1] + a1[1][1] * a2[1][1];
    return (t3);
}

/// The tensor multiplies itself by a coefficient
void SelfMultiply(tensor2D * t, double factor)
{
	double (*A)[2] = t->A;
    A[0][0] *= factor;
    A[0][1] *= factor;
    A[1][0] *= factor;
    A[1][1] *= factor;
}

/// Tensor divided by real number
tensor2D Tensor2DDivideCValue(tensor2D tensor, double factor)
{
    tensor2D result;
    double(*a1)[2] = tensor.A, (*a2)[2] = result.A;
    a2[0][0] = a1[0][0] / factor;
    a2[0][1] = a1[0][1] / factor;
    a2[1][0] = a1[1][0] / factor;
    a2[1][1] = a1[1][1] / factor;
    return (result);
}


/// Traces (i.e., the sum of the diagonals elements)
double Trace(tensor2D t)
{
    return (t.A[0][0] + t.A[1][1]);
}

/// Average of traces (i.e., the sum of the diagonals divided by the number of diagonal elements)
double TraceAverage(tensor2D t)
{
    return ((t.A[0][0] + t.A[1][1]) / 2.0);
}

/// Returns the maximum value of the tensor (i.e., max{ a[i][j] }, i,j=1,2)
double GetTensor2DMaxValue(tensor2D t)
{
    double maxval1 = max(t.A[0][0], t.A[0][1]);
    double maxval2 = max(t.A[1][0], t.A[1][1]);
    
    return( max(maxval1, maxval2) );
}

/// Convert C to Tensor
tensor2D CValueToTensor2D(double c)
{
	tensor2D tensor;
	double (*a)[2] = tensor.A;
	a[0][0] = c;
	a[0][1] = 0;
	a[1][0] = 0;
	a[1][1] = c;
	return(tensor);
}

/// Calculate strain rate
double CalcStrainRate(vec2D x1, vec2D x2, vec2D x3, vec2D v1, vec2D v2, vec2D v3)
{
    double strainRate = (Cross(Vec2DSub(v2,v1), Vec2DSub(x3, x1)) + Cross(Vec2DSub(x2, x1), Vec2DSub(v3, v1)))
                        / ( Cross(Vec2DSub(x2, x1), Vec2DSub(x3, x1)) );
	return(strainRate);
}

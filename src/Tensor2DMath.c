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


tensor2D DisplaceTensor(vec2D x2, vec2D x3, vec2D x2p, vec2D x3p)
{
    tensor2D tensor;
    double (*u)[2] = tensor.A;
    u[0][0] = (x3.Y * x2p.X - x2.Y * x3p.X) / (x3.Y * x2.X - x2.Y * x3.X);
    u[0][1] = (x3.X * x2p.X - x2.X * x3p.X) / (x3.X * x2.Y - x2.X * x3.Y);
    u[1][0] = (x3.Y * x2p.Y - x2.Y * x3p.Y) / (x3.Y * x2.X - x2.Y * x3.X);
    u[1][1] = (x3.X * x2p.Y - x2.X * x3p.Y) / (x3.X * x2.Y - x2.X * x3.Y);
    return (tensor);
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

/// Decompose tensor t1 into StrainRation + RotateRatio.
/// where StrainRation is a symmetric tensor,
///                   ┌     ┐
/// RotateRatio = vT *│0, -1│
///                   │1,  0│
///                   └     ┘
void TensorDepartion(tensor2D t, tensor2D * Strain, tensor2D * Rotate)
{
	double (*A)[2] = t.A;
    double vT = (A[1][0] - A[0][1]) / 2;
    Rotate->A[0][0] = 0;
    Rotate->A[0][1] = -vT;
    Rotate->A[1][0] = vT;
    Rotate->A[1][1] = 0;
    *Strain = Tensor2DSub(t, *Rotate);
}

/// Calculate the strain rate tensor from the two side vectors of the triangle and the velocity of motion
tensor2D CalcStrainRateTensor(vec2D x1, vec2D x2, vec2D v1, vec2D v2)
{
    tensor2D displaceTensor = DisplaceTensor(x1, x2, v1, v2);
    tensor2D strainRatio, rotateRatio;
    TensorDepartion(displaceTensor, &strainRatio, &rotateRatio);
    return (strainRatio);
}

/// Calculate strain rate
double CalcStrainRate(vec2D x1, vec2D x2, vec2D x3, vec2D v1, vec2D v2, vec2D v3)
{
    double strainRate = (Cross(Vec2DSub(v2,v1), Vec2DSub(x3, x1)) + Cross(Vec2DSub(x2, x1), Vec2DSub(v3, v1)))
        / (Cross(Vec2DSub(x2, x1), Vec2DSub(x3, x1)) + Cross(Vec2DSub(x2, x1), Vec2DSub(x3, x1)));
	return(strainRate);
}

/// Diagonalization of symmetric tensor
tensor2D Diagonalize(tensor2D t)
{
	double (*A)[2] = t.A;
    tensor2D diagTensor;
    /////////////////////////////////////////////////
    if (fabs(A[0][1]) < 1.0e-15 || fabs(A[0][1]) < 1.0e-10 * (fabs(A[0][0]) + fabs(A[1][1])))
    {
        diagTensor.A[0][0] = A[0][0];
        diagTensor.A[1][1] = A[1][1];
    }
    else
    {
        double b = (A[1][1] - A[0][0]) / A[0][1];
        double tanTheta = (-b + sqrt(b * b + 4.0)) / 2.0;
        double sinTheta = tanTheta / sqrt(1 + tanTheta * tanTheta);
        if (tanTheta < 1e-15)
        {
            diagTensor.A[0][0] = A[0][0];
            diagTensor.A[1][1] = A[1][1];
        }
        else
        {
            double cosTheta = sinTheta / tanTheta;
            /////////////////////////////////////////////////
            diagTensor.A[0][0] = A[0][0] * cosTheta * cosTheta + A[1][1] * sinTheta * sinTheta - 2 * A[0][1] * sinTheta * cosTheta;
            diagTensor.A[1][1] = A[0][0] * sinTheta * sinTheta + A[1][1] * cosTheta * cosTheta + 2 * A[0][1] * sinTheta * cosTheta;
        }
    }
    /////////////////////////////////////////////////
    return (diagTensor);
}
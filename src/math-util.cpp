
#include"math_util.hpp"

double dotproduct_Vec3(double a[], double b[])
{
    double product = 0.0;
    for(int i = 0; i < 3; ++i)
    {
        product += a[i] * b[i];
    }
    return product;
}

double sq(double a)
{
    return a*a;
}

double cb(double a)
{
    return a*a*a;
}

double v_sqr(double u, double v, double w)
{
    return u*u+v*v+w*w;
}
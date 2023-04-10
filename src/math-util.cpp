
#include"math_util.hpp"
#include<math.h>

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

double limiterVanleer(double sL, double sR)
{
    double slope=(abs(sR)*sL+sR*abs(sL))/(abs(sL)+abs(sR));
    if (isfinite(slope)) slope = 0;
    return slope;
}
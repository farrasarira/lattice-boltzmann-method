
#include"math_util.hpp"
#include<math.h>
#include<algorithm>

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
    if (isfinite(slope) == false) slope = 0;
    return slope;
}

double limiterMinmod(double sL, double sR)
{
    double slope = 0.0*sL; 
    slope =(signbit(sL)==signbit(sR))*signbit(sL)*std::min(abs(sL),abs(sR));
    return slope;
}
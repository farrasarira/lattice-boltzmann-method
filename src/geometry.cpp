#include "./headers/lbm.hpp"
#include <omp.h>

void cylinder_generator(LBM &lb, int D)
{
    int Nx = lb.getNx();
    int Ny = lb.getNy();
    int Nz = lb.getNz();
    double rad=D/2;         // radius of cylinder
    double cx=dx*(Nx/5);    // x-center of cylinder
    double cy=dx*(Ny/2);    // y-center of cylinder

    #pragma omp parallel for
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            for(int k=0; k<Nz; ++k)
            {
                double rx=dx*i-dx*0.5;
                double ry=dx*j-dx*0.5;
                if(sqrt(pow(cx-rx,2)+pow(cy-ry,2))<=rad){
                    lb.fluid1[i][j][k].type=TYPE_S;
                }
            }
        }
    }
};

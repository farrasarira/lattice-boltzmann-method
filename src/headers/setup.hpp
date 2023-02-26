#ifndef SETUP_H
    #define SETUP_H

    #include "lbm.hpp"
    #include "geometry.hpp"
    #include<iostream>

    LBM main_setup() // 2D Flow over cylinder --------------------------------------------------------
    {
        LBM lb(1000,500,1,0.002);
        int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

        int D = Ny/10; 
        double nu = 0.002;
        double Re = 1000;
        double u_max = Re * nu / D;
        double Ma = u_max * cs;

        cylinder_generator(lb, D);

        for(int i = 0; i < Nx ; ++i)
        {
            for(int j = 0; j < Ny; ++j)
            {
                for(int k = 0; k < Nz; ++k)
                {
                    if (lb.fluid1[i][j][k].type == TYPE_F)
                    {
                        lb.fluid1[i][j][k].rho = 1.0;
                        lb.fluid1[i][j][k].u = u_max;
                        lb.fluid1[i][j][k].v = 0;
                        lb.fluid1[i][j][k].w = 0;
                    }
                }
            }
        }

        return lb;
    }
    
#endif


#ifndef SETUP_H
    #define SETUP_H

    #include "lbm.hpp"
    #include "geometry.hpp"
    #include <iostream>
    #include <math.h>
    #include <omp.h>

    #if defined CYLINDER_2D
    LBM main_setup() // 2D Flow over cylinder --------------------------------------------------------
    {
        LBM lb(NX,NY,1,NU);
        int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

        int D = Ny/10; 
        double nu = lb.getNu();
        double u_max = RE * nu / D;

        cylinder_generator(lb, D);

        #pragma omp parallel for
        for(int i = 0; i < Nx ; ++i)
        {
            for(int j = 0; j < Ny; ++j)
            {
                for(int k = 0; k < Nz; ++k)
                {
                    if (j == 0 || j == Ny-1) 
                    {
                        lb.fluid1[i][j][k].type = TYPE_S;  // Solid Boundary in upper and lower side
                    }
                    if (i == 0 || i == Nx-1)
                    {
                        lb.fluid1[i][j][k].type = TYPE_E;  // Equilibrium inlet and outlet
                    }
                    if (k == 0 || k == Nz-1)
                    {
                        lb.fluid1[i][j][k].type = TYPE_P;
                    }
                    if (lb.fluid1[i][j][k].type == TYPE_F || lb.fluid1[i][j][k].type == TYPE_E)  // Fluid Domain
                    {
                        lb.fluid1[i][j][k].rho = 1.0;
                        lb.fluid1[i][j][k].u = u_max;
                        lb.fluid1[i][j][k].v = 0;
                        lb.fluid1[i][j][k].w = 0;
                    }
                }
            }
        }
        std::cout << "setup selesai" << std::endl;
        return lb;
    }

    #elif defined TAYLOR_GREEN_3D
    LBM main_setup() // 3D Taylor-Green Vortex
    {
        LBM lb(NX,NY,NZ,NU);
        int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

        const double u_max = RE * NU / (NX/(2.0*M_PI));    // maximum initial velocity
        std::cout << "umax      : " << u_max << std::endl;
        double periodicity = 1.0;
        const double a = (double)NX/periodicity;
        const double b = (double)NY/periodicity;
        const double c = (double)NZ/periodicity;
        
        #pragma omp parallel for
        for(int i = 0; i < Nx ; ++i)
        {
            for(int j = 0; j < Ny; ++j)
            {
                for(int k = 0; k < Nz; ++k)
                {
                    if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                    {
                        lb.fluid1[i][j][k].type = TYPE_P;
                    }
                    if (lb.fluid1[i][j][k].type == TYPE_F)
                    {
                        lb.fluid1[i][j][k].u =  u_max * sin(2.0*M_PI/a*(double)i) * cos(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                        lb.fluid1[i][j][k].v = -u_max * cos(2.0*M_PI/a*(double)i) * sin(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                        lb.fluid1[i][j][k].w = 0.0;
                        lb.fluid1[i][j][k].rho = 1.0 - u_max*u_max * 3.0/16.0 * (cos(4.0*M_PI/a*(double)i) + cos(4.0*M_PI/b*(double)j)) * (cos(4.0*M_PI/c*(double)k) + 2);
                        
                    }
                }
            }
        }
        return lb;
    }
    #elif defined TAYLOR_GREEN_2D
    LBM main_setup() // 2D Taylor-Green Vortex
    {
        LBM lb(NX,NY,1,NU);
        int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

        const double u_max = RE * NU / NX;    // maximum initial velocity
        std::cout << "umax      : " << u_max << std::endl;
        const unsigned int periodicity = 1;
        const double a = (double)NX/periodicity;
        const double b = (double)NX/periodicity;
        
        #pragma omp parallel for
        for(int i = 0; i < Nx ; ++i)
        {
            for(int j = 0; j < Ny; ++j)
            {
                for(int k = 0; k < Nz; ++k)
                {
                    if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1) // set periodic boundary condition
                    {
                        lb.fluid1[i][j][k].type = TYPE_P;
                    }
                    if (lb.fluid1[i][j][k].type == TYPE_F)
                    {
                        double X = i ;
                        double Y = j ;
                        lb.fluid1[i][j][k].u = -u_max * cos(2.0*M_PI/a*X) * sin(2.0*M_PI/b*Y) ;
                        lb.fluid1[i][j][k].v =  u_max * sin(2.0*M_PI/a*X) * cos(2.0*M_PI/b*Y) ;
                        lb.fluid1[i][j][k].rho = 1.0 - u_max*u_max * 3.0/4.0 * (cos(4.0*M_PI/a*X) + cos(4.0*M_PI/b*Y));
                        
                    }
                }
            }
        }
        return lb;
    }

    #elif defined CHANNEL_FLOW_3D
    LBM main_setup() // 3D Channel Flow
    {
        LBM lb(NX,NY,NY,NU);
        int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

        const double u_max = RE * NU / NY;    // maximum initial velocity
        std::cout << "umax      : " << u_max << std::endl;
        
        #pragma omp parallel for
        for(int i = 0; i < Nx ; ++i)
        {
            for(int j = 0; j < Ny; ++j)
            {
                for(int k = 0; k < Nz; ++k)
                {
                    if (j == 0 || j == Ny-1) 
                    {
                        lb.fluid1[i][j][k].type = TYPE_S;  // Solid Boundary in upper and lower side
                    }
                    if (i == 0 || i == Nx-1)
                    {
                        lb.fluid1[i][j][k].type = TYPE_E;  // Equilibrium inlet and outlet
                    }
                    if (k == 0 || k == Nz-1)
                    {
                        lb.fluid1[i][j][k].type = TYPE_P;
                    }
                    if (lb.fluid1[i][j][k].type == TYPE_F || lb.fluid1[i][j][k].type == TYPE_E)  // Fluid Domain
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

    #elif defined CYLINDER_3D
    LBM main_setup() // 3D Flow over cylinder --------------------------------------------------------
    {
        LBM lb(NX,NY,NZ,NU);
        int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

        int D = NY/10; 
        double nu = lb.getNu();
        double u_max = RE * nu / D;

        cylinder_generator(lb, D);

        #pragma omp parallel for
        for(int i = 0; i < Nx ; ++i)
        {
            for(int j = 0; j < Ny; ++j)
            {
                for(int k = 0; k < Nz; ++k)
                {
                    if (j == 0 || j == Ny-1 || k == 0 || k == Nz-1 ) 
                    {
                        lb.fluid1[i][j][k].type = TYPE_S;  // Solid Boundary in upper and lower side
                    }
                    if (i == 0 || i == Nx-1)
                    {
                        lb.fluid1[i][j][k].type = TYPE_E;  // Equilibrium inlet and outlet
                    }
                    if (lb.fluid1[i][j][k].type == TYPE_F || lb.fluid1[i][j][k].type == TYPE_E)  // Fluid Domain
                    {
                        lb.fluid1[i][j][k].rho = 1.0;
                        lb.fluid1[i][j][k].u = u_max;
                        lb.fluid1[i][j][k].v = 0;
                        lb.fluid1[i][j][k].w = 0;
                    }
                }
            }
        }
        std::cout << "setup selesai" << std::endl;
        return lb;
    }

    #endif
   
#endif


#include "lbm.hpp"
#include "geometry.hpp"
#include "units.hpp"
#include "cantera.hpp"
#include <iostream>
#include <math.h>
#include <omp.h>

Units_Conv units;

#if defined CYLINDER_2D
void main_setup() // 2D Flow over cylinder --------------------------------------------------------
{
    double RE = 1000;
    LBM lb(1000, 500, 1, 0.0025);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    double nu = lb.get_nu();

    int D = Ny/10; 
    cylinder_generator(lb, D);

    double u_max = RE*nu/D;    // [Velocity in Lattice Unit]
    double T0 = 0.03;       // [Temperature in Lattice Unit]

    // double si_u_max = 0.05; // [m/s]
    // double si_rho = 1.225; // [kg/m^3]
    // double si_temp = 273.15; // [K]
    // double si_cylinder_temp = 500.0; // [K]

    std::cout << "umax : " << u_max << std::endl;

    // auto sol = Cantera::newSolution("gri30.yaml","gri30");
    // auto gas = sol->thermo();
    // auto trans = sol->transport();
    // gas->setState_TR(si_temp, si_rho);
    
    // double si_D = RE*trans->viscosity()/si_rho/si_u_max; // [m]

    // units.set_m_kg_s(D, u_max, rho0, T0, si_D, si_u_max, si_rho, si_temp); // setting the conversion factor 

    // std::cout << gas->report() << std::endl;
    // std::cout << units.nu(trans->viscosity()/si_rho) << std::endl;

    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (j == 0 || j == Ny-1) 
                {
                    lb.mixture[i][j][k].type = TYPE_A;  // Solid Boundary in upper and lower side
                }
                if (i == 0 )
                {
                    lb.mixture[i][j][k].type = TYPE_I;  // Equilibrium inlet and outlet
                    lb.mixture[i][j][k].p = T0;
                    lb.mixture[i][j][k].temp = T0;
                    lb.mixture[i][j][k].u = u_max;
                }

                if (i == Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_O;
                    lb.mixture[i][j][k].p = T0;
                    lb.mixture[i][j][k].temp = T0;
                    
                }

                if (k == 0 || k == Nz-1)
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)  // Fluid Domain
                {
                    lb.mixture[i][j][k].p = T0;
                    lb.mixture[i][j][k].u = 0;
                    lb.mixture[i][j][k].v = 0;
                    lb.mixture[i][j][k].w = 0;
                    lb.mixture[i][j][k].temp = T0;
                }

            }
        }
    }

    lb.run(1000, 10);
}

#elif defined TAYLOR_GREEN_3D
void main_setup() // 3D Taylor-Green Vortex
{
    double RE = 10;
    LBM lb(100, 100, 100, 0.1);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    int NX = lb.get_NX(); int NY = lb.get_NY(); int NZ = lb.get_NZ();
    double NU = lb.get_nu();

    double periodicity = 1.0;
    const double a = (double)NX/periodicity;
    const double b = (double)NY/periodicity;
    const double c = (double)NZ/periodicity;

    double u_max = RE * NU / (NX/(2.*M_PI));
    std::cout << "umax : " << u_max << std::endl; 

    
    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].u =  u_max * sin(2.0*M_PI/a*(double)i) * cos(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].v = -u_max * cos(2.0*M_PI/a*(double)i) * sin(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].rho = 1.0 - u_max*u_max * 3.0/16.0 * (cos(4.0*M_PI/a*(double)i) + cos(4.0*M_PI/b*(double)j)) * (cos(4.0*M_PI/c*(double)k) + 2);
                    lb.mixture[i][j][k].temp = 0.1;
                }
            }
        }
    }
    
    lb.run(1000,10);
}
#elif defined TAYLOR_GREEN_2D
void main_setup() // 2D Taylor-Green Vortex
{
    double RE = 10;

    LBM lb(30, 30, 1, 0.01);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    double NX = lb.get_NX();

    double NU = lb.get_nu();

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
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    double X = i ;
                    double Y = j ;
                    lb.mixture[i][j][k].u = -u_max * cos(2.0*M_PI/a*X) * sin(2.0*M_PI/b*Y) ;
                    lb.mixture[i][j][k].v =  u_max * sin(2.0*M_PI/a*X) * cos(2.0*M_PI/b*Y) ;
                    lb.mixture[i][j][k].rho = 1.0 - u_max*u_max * 3.0/4.0 * (cos(4.0*M_PI/a*X) + cos(4.0*M_PI/b*Y));
                    
                }
            }
        }
    }
    
    lb.run(1000,10);
}

#elif defined CHANNEL_FLOW_3D
void main_setup() // 3D Channel Flow
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
                    lb.mixture[i][j][k].type = TYPE_S;  // Solid Boundary in upper and lower side
                }
                if (i == 0 || i == Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_E;  // Equilibrium inlet and outlet
                }
                if (k == 0 || k == Nz-1)
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F || lb.mixture[i][j][k].type == TYPE_E)  // Fluid Domain
                {
                    lb.mixture[i][j][k].rho = 1.0;
                    lb.mixture[i][j][k].u = u_max;
                    lb.mixture[i][j][k].v = 0;
                    lb.mixture[i][j][k].w = 0;
                }
            }
        }
    }
    return lb;
}

#elif defined CYLINDER_3D
void main_setup() // 3D Flow over cylinder --------------------------------------------------------
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
                    lb.mixture[i][j][k].type = TYPE_S;  // Solid Boundary in upper and lower side
                }
                if (i == 0 || i == Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_E;  // Equilibrium inlet and outlet
                }
                if (lb.mixture[i][j][k].type == TYPE_F || lb.mixture[i][j][k].type == TYPE_E)  // Fluid Domain
                {
                    lb.mixture[i][j][k].rho = 1.0;
                    lb.mixture[i][j][k].u = u_max;
                    lb.mixture[i][j][k].v = 0;
                    lb.mixture[i][j][k].w = 0;
                }
            }
        }
    }
    std::cout << "setup selesai" << std::endl;
    return lb;
}

#elif defined VISCOSITY_TEST
void main_setup() // 2D Viscos Test --------------------------------------------------------
{
    double NX = 1;
    double NY = 200;
    double NZ = 1;
    double nu = 0.1;

    LBM lb(NX,NY,NZ, nu);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    const double mach = 0.9;
    const double gamma = 1.4;
    const double a0 = 0.001; // amplitude
    const double T0 = 0.1;
    const double u0 = mach*sqrt(gamma*T0); 
    const double rho0 = 1.0;

    std::cout << "mach  : " << mach << std::endl;
    std::cout << "u0    : " << u0 << std::endl;

    lb.set_gamma(gamma);

    // double si_u_max = 1.0; // [m/s]
    // double si_rho = 1.225; // [kg/m^3]
    // double si_temp = 288.15; // [K]

    // auto gas = sols[0]->thermo();
    // auto trans = sols[0]->transport();
    // gas->setState_TR(si_temp, si_rho);

    // units.set_m_kg_s(NY, a0, rho0, T0, 1.0, si_u_max, si_rho, si_temp); // setting the conversion factor 

    // double si_pressure = gas->pressure();
    
    // double si_gas_constant = si_pressure/(si_rho*si_temp);
    // double gas_const = units.R(si_gas_constant);
    
    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    //double X = i ;
                    double Y = j ;
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0 ;
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].p = 0.1;
                    double rho = rho0 + a0*sin(2.0*M_PI/NY*Y);
                    lb.mixture[i][j][k].temp = lb.mixture[i][j][k].p / rho;
                    
                }
            }
        }
    }
    std::cout << "u : " << lb.mixture[Nx/2][Ny/2][Nz/2].u << std::endl;
    lb.run(10000,1000);
}

#elif defined SOD_SHOCK
void main_setup() // Sos shock tube test case --------------------------------------------------------
{
    int NX = 3000; 
    int NY = 1; 
    int NZ = 1;
    
    double NU = 0.025;

    double si_len = 1E-0;    // [m]
    double si_u_max = 10.0;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 300.0;// [K]

    units.set_m_kg_s(NX, VEL0, RHO0, TEMP0, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    LBM lb(NX, NY, NZ, NU);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    
    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (i==0 || i==Nx-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {   
                    lb.mixture[i][j][k].p = smooth(2.0, 0.5, i, 0.5*Nx, 0.3);   
                    // std::cout << " p : " << units.si_p(2.0) << std::endl;
                    lb.mixture[i][j][k].temp = smooth(0.1, 0.1, i, 0.5*Nx, 0.3);   
                    if ((float)i/(float)Nx <= 0.5 )
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                        // lb.mixture[i][j][k].p = 0.01;
                        // lb.mixture[i][j][k].temp = 0.2;
                    }
                    else
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                        // lb.mixture[i][j][k].p = 0.001;
                        // lb.mixture[i][j][k].temp = 0.025;
                    }                       
                }
            }
        }
    }

    lb.run(1000,100);
}


#elif defined SOD_SHOCK_SIUNIT
void main_setup() // SOD SHOCK TUBE WITH SI UNIT --------------------------------------------------------
{
    std::cout << "-------------------- SOD SCHOCK SIUNIT ---------------------" << std::endl;
    int NX = 3000; 
    int NY = 1; 
    int NZ = 1;
    
    double si_len = 1E-2;    // [m]
    double si_u_max = 1E+3;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 300.0;// [K]

    units.set_m_kg_s(NX, VEL0, RHO0, 0.1, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    std::vector<std::string> species = { "Ar" };
        
    LBM lb(NX, NY, NZ, species);
    lb.set_diffusionModel("Stefan-Maxwell");
    // lb.set_diffusionModel("Mixture-Averaged");
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (i==0 || i==Nx-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.species[0][i][j][k].X = 1.0;

                    // lb.mixture[i][j][k].p = smooth(2*units.p(Cantera::OneAtm), units.p(Cantera::OneAtm), i, 0.5*Nx, 0.3);      
                    // lb.mixture[i][j][k].temp = smooth(units.temp(400.0), units.temp(400.0), i, 0.5*Nx, 0.3);                  

                    lb.mixture[i][j][k].p = smooth(2*units.p(Cantera::OneAtm), units.p(Cantera::OneAtm), i, 0.5*Nx, 0.3);      
                    lb.mixture[i][j][k].temp = smooth(units.temp(300.0), units.temp(300.0), i, 0.5*Nx, 0.3); 

                    if ((float)i/(float)Nx <= 0.5 )
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0; 
                    }
                    else 
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                    }                       
                }
            }
        }
    }

    lb.run(10000,1000);
}

#elif defined COUETTE_FLOW

void main_setup() // 2D Couette Flow --------------------------------------------------------
{
    int NX = 2; 
    int NY = 100; 
    int NZ = 1;
    
    double si_len = 1E-3;    // [m]
    double si_u_max = 1.0E+2;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 300.0;// [K]
    double nu = 0.025;

    units.set_m_kg_s(NX, VEL0, RHO0, TEMP0, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    
    LBM lb(NX, NY, NZ, nu);

    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();


    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( i==0 || i==Nx-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (j==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].temp = 0.1;
                }
                if (j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = 0.05;
                    lb.mixture[i][j][k].temp = 0.1003125;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].temp = 0.1;
                    lb.mixture[i][j][k].p = 1.0;
                    lb.mixture[i][j][k].u = 0.0;
                }
            }
        }
    }

    lb.set_prtl(0.25);
    lb.set_gasconst(0.4/1.4);
    lb.run(200000,10000);
}

#endif
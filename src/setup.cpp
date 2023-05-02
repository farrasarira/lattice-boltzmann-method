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
    LBM lb(NX,NY,1);
    int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

    int D = Ny/10; 
    cylinder_generator(lb, D);

    double u_max = 0.01; // [lattice speed unit]
    double si_u_max = 0.05; // [m/s]
    double si_rho = 1.225; // [kg/m^3]
    double si_temp = 273.15; // [K]
    double si_cylinder_temp = 500.0; // [K]

    std::cout << "umax : " << u_max << std::endl;

    sol = Cantera::newSolution("gri30.yaml","gri30");
    gas = sol->thermo();
    trans = sol->transport();
    gas->setState_TR(si_temp, si_rho);
    
    double si_D = RE*trans->viscosity()/si_rho/si_u_max; // [m]

    units.set_m_kg_s(D, u_max, rho0, TREF, si_D, si_u_max, si_rho, si_temp); // setting the conversion factor 

    std::cout << gas->report() << std::endl;
    std::cout << units.nu(trans->viscosity()/si_rho) << std::endl;

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
                    lb.fluid1[i][j][k].temp = TREF;
                }
                if(lb.fluid1[i][j][k].type == TYPE_S)
                {
                    lb.fluid1[i][j][k].temp = units.temp(si_cylinder_temp);
                }
            }
        }
    }
    return lb;
}

#elif defined TAYLOR_GREEN_3D
void main_setup() // 3D Taylor-Green Vortex
{
    LBM lb(NX,NY,NZ);
    int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

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
                    lb.fluid1[i][j][k].type = TYPE_P;
                }
                if (lb.fluid1[i][j][k].type == TYPE_F)
                {
                    lb.fluid1[i][j][k].u =  u_max * sin(2.0*M_PI/a*(double)i) * cos(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.fluid1[i][j][k].v = -u_max * cos(2.0*M_PI/a*(double)i) * sin(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.fluid1[i][j][k].w = 0.0;
                    lb.fluid1[i][j][k].rho = 1.0 - u_max*u_max * 3.0/16.0 * (cos(4.0*M_PI/a*(double)i) + cos(4.0*M_PI/b*(double)j)) * (cos(4.0*M_PI/c*(double)k) + 2);
                    lb.fluid1[i][j][k].temp = 0.1;
                }
            }
        }
    }
    return lb;
}
#elif defined TAYLOR_GREEN_2D
void main_setup() // 2D Taylor-Green Vortex
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

#elif defined VISCOSITY_TEST
void main_setup() // 2D Viscos Test --------------------------------------------------------
{
    LBM lb(NX,NY,NZ);
    int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

    const double mach = 0.9;
    const double gamma = GAMMA;
    const double a0 = 0.001; // amplitude
    const double T0 = 0.1;
    const double u0 = mach*sqrt(gamma*T0); 
    const double rho0 = 1.0;

    std::cout << "mach  : " << mach << std::endl;
    std::cout << "u0    : " << u0 << std::endl;

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
                    lb.fluid1[i][j][k].type = TYPE_P;
                }
                if (lb.fluid1[i][j][k].type == TYPE_F)
                {
                    //double X = i ;
                    double Y = j ;
                    lb.fluid1[i][j][k].u = a0*sin(2.0*M_PI/NY*Y);
                    lb.fluid1[i][j][k].v = u0 ;
                    lb.fluid1[i][j][k].w = 0.0;
                    lb.fluid1[i][j][k].rho = rho0;
                    lb.fluid1[i][j][k].temp = T0;
                    
                }
            }
        }
    }
    std::cout << "u : " << lb.fluid1[Nx/2][Ny/2][Nz/2].u << std::endl;
    return lb;
}

#elif defined SOD_SHOCK
void main_setup() // 2D Viscos Test --------------------------------------------------------
{
    int NX = 3000; 
    int NY = 5; 
    int NZ = 5;
    double NU = 0.025;
    
    double si_len = 0.1;    // [m]
    double si_u_max = 1.0;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 288.15;// [K]

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
                    lb.mixture[i][j][k].type = TYPE_S;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    if ((float)i/(float)Nx <= 0.5 )
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                        lb.mixture[i][j][k].rho = units.rho(0.7);
                        lb.mixture[i][j][k].temp = units.temp(300.0);
                    }
                    else
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                        lb.mixture[i][j][k].rho = units.rho(1.2);
                        lb.mixture[i][j][k].temp = units.temp(300.0);
                    }                       
                }
            }
        }
    }

    std::cout << "units.temp : " << units.temp(300.0) << std::endl;
    lb.run(1000,10);
}

#elif defined TERNARY_DIFFUSION

void main_setup() // 2D Viscos Test --------------------------------------------------------
{
    int NX = 900; 
    int NY = 5; 
    int NZ = 5;
    
    double si_len = 0.001;    // [m]
    double si_u_max = 1000.0;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 300.0;// [K]

    units.set_m_kg_s(NX, VEL0, RHO0, TEMP0, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    std::vector<std::string> species = { "H2", "AR", "CH4" };
    
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz(); int nSpecies = lb.get_nSpecies();

    auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();
    std::vector<double> XL (gas->nSpecies());
    std::vector<double> XR (gas->nSpecies()); 

    XL[gas->speciesIndex("H2")]  = 0.491;
    XL[gas->speciesIndex("AR")]  = 0.509;
    XL[gas->speciesIndex("CH4")] = 0.0;

    XR[gas->speciesIndex("H2")]  = 0.0;
    XR[gas->speciesIndex("AR")]  = 0.485;
    XR[gas->speciesIndex("CH4")] = 0.515;

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
                    lb.mixture[i][j][k].type = TYPE_S;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm);
                    lb.mixture[i][j][k].temp = units.temp(300.0);
                    if ((float)i/(float)Nx <= 0.5 )
                    {
                        for (int a = 0; a < nSpecies; ++a) lb.species[a][i][j][k].X = XL[gas->speciesIndex(species[a])];
                    }
                    else
                    {
                        for (int a = 0; a < nSpecies; ++a) lb.species[a][i][j][k].X = XR[gas->speciesIndex(species[a])];
                    }        
                }
            }
        }
    }


    lb.run(1000,10);
}

#endif



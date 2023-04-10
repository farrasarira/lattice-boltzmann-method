
#include "cantera.hpp"
#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include <omp.h>
#include <numeric>
#include <vector>

// using namespace Cantera;

// std::vector<std::shared_ptr<Cantera::Solution>> sols;

LBM::LBM(int Nx, int Ny, int Nz)
{
    this->Nx = Nx + 2;
    this->Ny = Ny + 2;
    this->Nz = Nz + 2;
    // allocate memory for lattice
    fluid1 = new LATTICE **[this->Nx];
    for (int i = 0; i < this->Nx; ++i)
    {
        fluid1[i] = new LATTICE *[this->Ny];
        for (int j = 0; j < this->Ny; ++j)
        {
            fluid1[i][j] = new LATTICE [this->Nz];
        }
    }

    // int nThreads  = omp_get_max_threads();

    // for(int i = 0; i < nThreads; ++i)
    // {
    //     auto sol = newSolution("gri30.yaml", "gri30");
    //     sols.emplace_back(sol);
    // }
}


void LBM::Init()
{
    #ifdef LBM_EXTEND
    std::vector<std::vector<std::vector<double>>> dQdevx(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevy(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevz(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevx_lim(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevy_lim(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevz_lim(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));

    #pragma omp parallel for //schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (fluid1[i][j][k].type == TYPE_F || fluid1[i][j][k].type == TYPE_E)     
                {
                    if (i == 1)
                    {
                        dQdevx[i][j][k] = ((fluid1[i+1][j][k].rho*fluid1[i+1][j][k].u*(1-3*fluid1[i+1][j][k].gas_const*fluid1[i+1][j][k].temp)-fluid1[i+1][j][k].rho*cb(fluid1[i+1][j][k].u)) - (fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].u)) ) /(dx);
                    }
                    else
                    {
                        dQdevx[i][j][k] = ((fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].u)) - (fluid1[i-1][j][k].rho*fluid1[i-1][j][k].u*(1-3*fluid1[i-1][j][k].gas_const*fluid1[i-1][j][k].temp)-fluid1[i-1][j][k].rho*cb(fluid1[i-1][j][k].u)) ) /(dx);
                    }

                    if (j == 1)
                    {
                        dQdevy[i][j][k] = ((fluid1[i][j+1][k].rho*fluid1[i][j+1][k].v*(1-3*fluid1[i][j+1][k].gas_const*fluid1[i][j+1][k].temp)-fluid1[i][j+1][k].rho*cb(fluid1[i][j+1][k].v)) - (fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].v))) /(dy);
                    }
                    else
                    {
                        dQdevy[i][j][k] = ((fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].v)) - (fluid1[i][j-1][k].rho*fluid1[i][j-1][k].v*(1-3*fluid1[i][j-1][k].gas_const*fluid1[i][j-1][k].temp)-fluid1[i][j-1][k].rho*cb(fluid1[i][j-1][k].v)) ) /(dy);
                    }

                    if (k == 1)
                    {
                        dQdevz[i][j][k] = ((fluid1[i][j][k+1].rho*fluid1[i][j][k+1].w*(1-3*fluid1[i][j][k+1].gas_const*fluid1[i][j][k+1].temp)-fluid1[i][j][k+1].rho*cb(fluid1[i][j][k+1].w)) - (fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].w)))  /(dz);
                    }
                    else
                    {
                        dQdevz[i][j][k] = ((fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].w)) - (fluid1[i][j][k-1].rho*fluid1[i][j][k-1].w*(1-3*fluid1[i][j][k-1].gas_const*fluid1[i][j][k-1].temp)-fluid1[i][j][k-1].rho*cb(fluid1[i][j][k-1].w)) ) /(dz);
                    }
                }
            }
        }
    }

    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (fluid1[i][j][k].type == TYPE_F || fluid1[i][j][k].type == TYPE_E)     
                {
                    if (i == Nx-2) dQdevx_lim[i][j][k] = limiterVanleer(dQdevx[i-1][j][k],dQdevx[i][j][k]);
                    else dQdevx_lim[i][j][k] = limiterVanleer(dQdevx[i][j][k],dQdevx[i+1][j][k]);

                    if (j == Ny-2) dQdevy_lim[i][j][k] = limiterVanleer(dQdevy[i][j-1][k],dQdevy[i][j][k]);
                    else dQdevy_lim[i][j][k] = limiterVanleer(dQdevy[i][j][k],dQdevy[i][j+1][k]);

                    if (k == Nz-2) dQdevz_lim[i][j][k] = limiterVanleer(dQdevz[i][j][k-1],dQdevz[i][j][k]);
                    else dQdevz_lim[i][j][k] = limiterVanleer(dQdevz[i][j][k],dQdevz[i][j][k+1]);
                }
            }
        }
    }


    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (fluid1[i][j][k].type == TYPE_F || fluid1[i][j][k].type == TYPE_E)     
                {   
                    // size_t rank = omp_get_thread_num();
                    // auto gas = sols[rank]->thermo();
                    // auto trans = sols[rank]->transport();
                    // gas->setState_TR(units.si_temp(fluid1[i][j][k].temp), units.si_rho(fluid1[i][j][k].rho));

                    // fluid1[i][j][k].p = fluid1[i][j][k].rho * fluid1[i][j][k].gas_const * fluid1[i][j][k].temp;

                    // fluid1[i][j][k].mu = units.mu(trans->viscosity());
                    // //fluid1[i][j][k].nu = NU; // debugg...
                    // fluid1[i][j][k].omega = 2*fluid1[i][j][k].p*dt_sim / (fluid1[i][j][k].p*dt_sim + 2*fluid1[i][j][k].mu);
        
                    // fluid1[i][j][k].conduc_coeff = units.thermalConductivity(trans->thermalConductivity());
                    // fluid1[i][j][k].cp= units.cp(gas->cp_mass());
                    // fluid1[i][j][k].omega1 = 2*fluid1[i][j][k].p*fluid1[i][j][k].cp*dt_sim / (fluid1[i][j][k].p*fluid1[i][j][k].cp*dt_sim + 2*fluid1[i][j][k].conduc_coeff);
                    // //double konst = PR * fluid1[i][j][k].omega / (2 - fluid1[i][j][k].omega); // debugg...
                    // //fluid1[i][j][k].omega1 = 2 * konst /(1 + konst); // debugg...

                    // fluid1[i][j][k].internalEnergy = units.energy_mass(gas->intEnergy_mass());
                    // //fluid1[i][j][k].internalEnergy = 3.0/2.0 * fluid1[i][j][k].gas_const * fluid1[i][j][k].temp; // debugg...
                    // fluid1[i][j][k].enthalpy = units.energy_mass(gas->enthalpy_mass());
                    // //fluid1[i][j][k].enthalpy = 5.0/2.0 * fluid1[i][j][k].gas_const * fluid1[i][j][k].temp; // debugg...
                    // fluid1[i][j][k].totalEnergy = fluid1[i][j][k].internalEnergy + v_sqr(fluid1[i][j][k].u, fluid1[i][j][k].v, fluid1[i][j][k].w)/2.0; 
                                          
                    fluid1[i][j][k].p = fluid1[i][j][k].rho * fluid1[i][j][k].gas_const * fluid1[i][j][k].temp;
                    fluid1[i][j][k].mu = NU*fluid1[i][j][k].rho;
                    fluid1[i][j][k].omega = 2*fluid1[i][j][k].p*dt_sim / (fluid1[i][j][k].p*dt_sim + 2*fluid1[i][j][k].mu);
                    
                    double cv = fluid1[i][j][k].gas_const / (GAMMA - 1.0);
                    double cp = cv + fluid1[i][j][k].gas_const;
                    fluid1[i][j][k].internalEnergy = cv * fluid1[i][j][k].temp;
                    fluid1[i][j][k].totalEnergy = fluid1[i][j][k].internalEnergy + 0.5 * v_sqr(fluid1[i][j][k].u, fluid1[i][j][k].v, fluid1[i][j][k].w);
                    fluid1[i][j][k].enthalpy = cp * fluid1[i][j][k].temp; // H = Cp * T = (Cv + 1) * T
                    fluid1[i][j][k].conduc_coeff = fluid1[i][j][k].mu*cp/PR;
                    fluid1[i][j][k].omega1 = 2*fluid1[i][j][k].p*cp*dt_sim / (fluid1[i][j][k].p*cp*dt_sim + 2*fluid1[i][j][k].conduc_coeff);

                    for (int l = 0; l < npop; ++l)
                    {
                        // ------------- Mass and Momentum Initialization -----------------------------
                        double P = 0.0;
                        double eps = 0.0;
                        double feq = fluid1[i][j][k].rho;

                            eps = fluid1[i][j][k].u;
                            P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + sq(fluid1[i][j][k].u) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*dQdevx_lim[i][j][k];
                            if (cx[l] == 0) feq *= (1 - P);
                            else if (cx[l] == 1) feq *= (eps+P)/2;
                            else if (cx[l] == -1) feq*= (-eps+P)/2;
                        
                        #if NDIM == 2 || NDIM == 3
                            eps = fluid1[i][j][k].v;
                            P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + sq(fluid1[i][j][k].v)+ dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*dQdevy_lim[i][j][k];
                            if (cy[l] == 0) feq *= (1 - P);
                            else if (cy[l] == 1) feq *= (eps+P)/2;
                            else if (cy[l] == -1) feq*= (-eps+P)/2;
                        #endif

                        #if NDIM == 3
                            eps = fluid1[i][j][k].w;
                            P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + sq(fluid1[i][j][k].w) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*dQdevz_lim[i][j][k];
                            if (cz[l] == 0) feq *= (1 - P);
                            else if (cz[l] == 1) feq *= (eps+P)/2;
                            else if (cz[l] == -1) feq*= (-eps+P)/2;
                        #endif
                    
                        fluid1[i][j][k].fpc[l]=feq;
                        fluid1[i][j][k].f[l]=feq;
                    }

                    // ------------------- Energy Initialization -----------------------------------
                    double theta =  fluid1[i][j][k].gas_const*fluid1[i][j][k].temp; //1./3.;//
                    double velocity[3] = {fluid1[i][j][k].u, fluid1[i][j][k].v, fluid1[i][j][k].w};
                    double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                {0., 0., 0.},
                                                {0., 0., 0.}};
                    // double eq_p_tensor_corr[3][3] = {{0., 0., 0.},    // pressure tensor
                    //                             {0., 0., 0.},
                    //                             {0., 0., 0.}};
                    memset(fluid1[i][j][k].p_tensor, 0.0, sizeof(fluid1[i][j][k].p_tensor));  

                    for (int l = 0; l < npop; ++l)
                    {
                        double velocity_set[3] = {cx[l], cy[l], cz[l]};
                        for(int p=0; p < 3; ++p)
                        {
                            for(int q=0; q < 3; ++q)
                            {
                                fluid1[i][j][k].p_tensor[p][q] += fluid1[i][j][k].f[l]*velocity_set[p]*velocity_set[q]; 
                                eq_p_tensor[p][q] = (p==q) ? fluid1[i][j][k].p+fluid1[i][j][k].rho*velocity[p]*velocity[q] : fluid1[i][j][k].rho*velocity[p]*velocity[q]; 
                                //eq_p_tensor_corr[p][q] += feq_temp[l]*velocity_set[p]*velocity_set[q];//
                            }  
                        }
                    }
                                       
                    for (int l = 0; l < npop; ++l)
                    {
                        double velocity_set[3] = {cx[l], cy[l], cz[l]};
                        double total_enthalpy = fluid1[i][j][k].enthalpy + 0.5 * v_sqr(fluid1[i][j][k].u,fluid1[i][j][k].v,fluid1[i][j][k].w);
                        double eq_heat_flux[3] = {total_enthalpy*fluid1[i][j][k].rho*fluid1[i][j][k].u,
                                                  total_enthalpy*fluid1[i][j][k].rho*fluid1[i][j][k].v,
                                                  total_enthalpy*fluid1[i][j][k].rho*fluid1[i][j][k].w};
                        double eq_R_tensor[3][3] = {{0., 0., 0.},    // second-order moment of g
                                                    {0., 0., 0.},
                                                    {0., 0., 0.}};
                        for(int p=0; p < 3; ++p)
                        {
                            for(int q=0; q < 3; ++q)
                            {
                                eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + fluid1[i][j][k].p*velocity[p]*velocity[q];
                            }
                        }
                        

                        double geq = fluid1[i][j][k].rho * fluid1[i][j][k].totalEnergy;

                        geq += dotproduct_Vec3(eq_heat_flux, velocity_set) / theta;

                        double matA[3][3] = {{eq_R_tensor[0][0]-fluid1[i][j][k].rho*fluid1[i][j][k].totalEnergy*theta, eq_R_tensor[0][1]                                                      , eq_R_tensor[0][2]},
                                             {eq_R_tensor[1][0]                                                      , eq_R_tensor[1][1]-fluid1[i][j][k].rho*fluid1[i][j][k].totalEnergy*theta, eq_R_tensor[1][2]},
                                             {eq_R_tensor[2][0]                                                      , eq_R_tensor[2][1]                                                      , eq_R_tensor[2][2]-fluid1[i][j][k].rho*fluid1[i][j][k].totalEnergy*theta}};
                        double matB[3][3] = {{cx[l]*cx[l]-theta  , cx[l]*cy[l]          , cx[l]*cz[l]},
                                             {cy[l]*cx[l]        , cy[l]*cy[l]-theta    , cy[l]*cz[l]},
                                             {cz[l]*cx[l]        , cz[l]*cy[l]          , cz[l]*cz[l]-theta}};
                        double result_AB = (matA[0][0]*matB[0][0] + matA[1][0]*matB[1][0] + matA[2][0]*matB[2][0]) + (matA[0][1]*matB[0][1] + matA[1][1]*matB[1][1] + matA[2][1]*matB[2][1]) + (matA[0][2]*matB[0][2] + matA[1][2]*matB[1][2] + matA[2][2]*matB[2][2]);
                        
                        geq += result_AB/(2.0 * theta * theta);
                        
                        double weight = 1.0;
                        for (int m = 0; m < 3; ++m)
                        {
                            if (velocity_set[m] == 0) weight *= (1 - theta);
                            else weight *= theta / 2.0;
                        }
                        geq = weight * geq;

                        double B;
                        double Z;
                        for (int m = 0; m < 3; ++m)
                        {
                            if (v_sqr(velocity_set[0], velocity_set[1], velocity_set[2]) == 0) B = 1;
                            else if (v_sqr(velocity_set[0], velocity_set[1], velocity_set[2]) == 1) B = -0.5*abs(velocity_set[m]);
                            else B = 0;
                            
                            Z = (1-3*theta)/(2*theta) * (eq_R_tensor[m][m]-theta*fluid1[i][j][k].rho*fluid1[i][j][k].totalEnergy);
                            
                            geq += B*Z; 
                        }

                        fluid1[i][j][k].g[l] = geq;
                        fluid1[i][j][k].gpc[l] = geq;    
                    }
                }
            }
        }
    }
   
    #else
    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (fluid1[i][j][k].type == TYPE_F || fluid1[i][j][k].type == TYPE_E)     
                {   
                    double uu = fluid1[i][j][k].u*fluid1[i][j][k].u + fluid1[i][j][k].v*fluid1[i][j][k].v + fluid1[i][j][k].w*fluid1[i][j][k].w;

                    for (int l = 0; l < npop; ++l)
                    {
                        double cu = cx[l]*fluid1[i][j][k].u + cy[l]*fluid1[i][j][k].v + cz[l]*fluid1[i][j][k].w;
                        double feq = wi[l]*fluid1[i][j][k].rho*(1.0+3.0*cu+4.5*cu*cu-1.5*uu);
                        fluid1[i][j][k].fpc[l]=feq;
                        fluid1[i][j][k].f[l]=feq;
                    }
                }
            }
        }
    }
    #endif
}

void LBM::Collide_BGK()
{
    std::vector<std::vector<std::vector<double>>> dQdevx(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevy(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevz(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevx_lim(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevy_lim(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevz_lim(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));

    #pragma omp parallel for //schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (fluid1[i][j][k].type == TYPE_F || fluid1[i][j][k].type == TYPE_E)     
                {
                    if (i == 1)
                    {
                        dQdevx[i][j][k] = ((fluid1[i+1][j][k].rho*fluid1[i+1][j][k].u*(1-3*fluid1[i+1][j][k].gas_const*fluid1[i+1][j][k].temp)-fluid1[i+1][j][k].rho*cb(fluid1[i+1][j][k].u)) - (fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].u)) ) /(dx);
                    }
                    else
                    {
                        dQdevx[i][j][k] = ((fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].u)) - (fluid1[i-1][j][k].rho*fluid1[i-1][j][k].u*(1-3*fluid1[i-1][j][k].gas_const*fluid1[i-1][j][k].temp)-fluid1[i-1][j][k].rho*cb(fluid1[i-1][j][k].u)) ) /(dx);
                    }

                    if (j == 1)
                    {
                        dQdevy[i][j][k] = ((fluid1[i][j+1][k].rho*fluid1[i][j+1][k].v*(1-3*fluid1[i][j+1][k].gas_const*fluid1[i][j+1][k].temp)-fluid1[i][j+1][k].rho*cb(fluid1[i][j+1][k].v)) - (fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].v))) /(dy);
                    }
                    else
                    {
                        dQdevy[i][j][k] = ((fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].v)) - (fluid1[i][j-1][k].rho*fluid1[i][j-1][k].v*(1-3*fluid1[i][j-1][k].gas_const*fluid1[i][j-1][k].temp)-fluid1[i][j-1][k].rho*cb(fluid1[i][j-1][k].v)) ) /(dy);
                    }

                    if (k == 1)
                    {
                        dQdevz[i][j][k] = ((fluid1[i][j][k+1].rho*fluid1[i][j][k+1].w*(1-3*fluid1[i][j][k+1].gas_const*fluid1[i][j][k+1].temp)-fluid1[i][j][k+1].rho*cb(fluid1[i][j][k+1].w)) - (fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].w)))  /(dz);
                    }
                    else
                    {
                        dQdevz[i][j][k] = ((fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*cb(fluid1[i][j][k].w)) - (fluid1[i][j][k-1].rho*fluid1[i][j][k-1].w*(1-3*fluid1[i][j][k-1].gas_const*fluid1[i][j][k-1].temp)-fluid1[i][j][k-1].rho*cb(fluid1[i][j][k-1].w)) ) /(dz);
                    }
                }
            }
        }
    }

    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (fluid1[i][j][k].type == TYPE_F || fluid1[i][j][k].type == TYPE_E)     
                {
                    if (i == Nx-2) dQdevx_lim[i][j][k] = limiterVanleer(dQdevx[i-1][j][k],dQdevx[i][j][k]);
                    else dQdevx_lim[i][j][k] = limiterVanleer(dQdevx[i][j][k],dQdevx[i+1][j][k]);

                    if (j == Ny-2) dQdevy_lim[i][j][k] = limiterVanleer(dQdevy[i][j-1][k],dQdevy[i][j][k]);
                    else dQdevy_lim[i][j][k] = limiterVanleer(dQdevy[i][j][k],dQdevy[i][j+1][k]);

                    if (k == Nz-2) dQdevz_lim[i][j][k] = limiterVanleer(dQdevz[i][j][k-1],dQdevz[i][j][k]);
                    else dQdevz_lim[i][j][k] = limiterVanleer(dQdevz[i][j][k],dQdevz[i][j][k+1]);
                }
            }
        }
    }

    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int k = 0; k < Nz; ++k)
            {
                if (fluid1[i][j][k].type==TYPE_F)
                {                   
                    #ifdef LBM_EXTEND

                    for (int l = 0; l < npop; ++l)
                    {
                        double P = 0.0;
                        double eps = 0.0;
                        double feq = fluid1[i][j][k].rho;

                            eps = fluid1[i][j][k].u;
                            P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + sq(fluid1[i][j][k].u) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*dQdevx_lim[i][j][k];
                            if (cx[l] == 0) feq *= (1 - P);
                            else if (cx[l] == 1) feq *= (eps+P)/2;
                            else if (cx[l] == -1) feq*= (-eps+P)/2;
                        
                        #if NDIM == 2 || NDIM == 3
                            eps = fluid1[i][j][k].v;
                            P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + sq(fluid1[i][j][k].v)+ dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*dQdevy_lim[i][j][k];
                            if (cy[l] == 0) feq *= (1 - P);
                            else if (cy[l] == 1) feq *= (eps+P)/2;
                            else if (cy[l] == -1) feq*= (-eps+P)/2;
                        #endif

                        #if NDIM == 3
                            eps = fluid1[i][j][k].w;
                            P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + sq(fluid1[i][j][k].w) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*dQdevz_lim[i][j][k];
                            if (cz[l] == 0) feq *= (1 - P);
                            else if (cz[l] == 1) feq *= (eps+P)/2;
                            else if (cz[l] == -1) feq*= (-eps+P)/2;
                        #endif

                        fluid1[i][j][k].fpc[l]=(1.0-fluid1[i][j][k].omega)*fluid1[i][j][k].f[l] + fluid1[i][j][k].omega*feq;

                    }

                    // ------------------- Energy Initialization -----------------------------------
                    double theta = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp;
                    double velocity[3] = {fluid1[i][j][k].u, fluid1[i][j][k].v, fluid1[i][j][k].w};
                    double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                {0., 0., 0.},
                                                {0., 0., 0.}};
                    memset(fluid1[i][j][k].p_tensor, 0.0, sizeof(fluid1[i][j][k].p_tensor));  
                    
                    for (int l = 0; l < npop; ++l)
                    {
                        double velocity_set[3] = {cx[l], cy[l], cz[l]};
                        for(int p=0; p < 3; ++p)
                        {
                            for(int q=0; q < 3; ++q)
                            {
                                fluid1[i][j][k].p_tensor[p][q] += fluid1[i][j][k].f[l]*velocity_set[p]*velocity_set[q]; 
                                eq_p_tensor[p][q] = ((p==q)? fluid1[i][j][k].p : 0) + fluid1[i][j][k].rho*velocity[p]*velocity[q];  
                            }  
                        }
                    }
                    
                    for (int l = 0; l < npop; ++l)
                    {
                        // ------------------- Energy Collision -----------------------------------
                        double velocity_set[3] = {cx[l], cy[l], cz[l]};
                        double total_enthalpy = fluid1[i][j][k].enthalpy+v_sqr(fluid1[i][j][k].u,fluid1[i][j][k].v,fluid1[i][j][k].w)/2.0;
                        double eq_heat_flux[3] = {total_enthalpy*fluid1[i][j][k].rho*fluid1[i][j][k].u,
                                                  total_enthalpy*fluid1[i][j][k].rho*fluid1[i][j][k].v,
                                                  total_enthalpy*fluid1[i][j][k].rho*fluid1[i][j][k].w};
                        
                        double str_heat_flux[3] ={fluid1[i][j][k].energy_flux[0] - velocity[0]*(fluid1[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) - velocity[1]*(fluid1[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) - velocity[2]*(fluid1[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) - 0.5*dt_sim*fluid1[i][j][k].u*dQdevx_lim[i][j][k],
                                                  fluid1[i][j][k].energy_flux[1] - velocity[0]*(fluid1[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) - velocity[1]*(fluid1[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) - velocity[2]*(fluid1[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) - 0.5*dt_sim*fluid1[i][j][k].v*dQdevy_lim[i][j][k],
                                                  fluid1[i][j][k].energy_flux[2] - velocity[0]*(fluid1[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) - velocity[1]*(fluid1[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) - velocity[2]*(fluid1[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) - 0.5*dt_sim*fluid1[i][j][k].w*dQdevz_lim[i][j][k]};
                        
                        double eq_R_tensor[3][3] = {{0., 0., 0.},    // second-order moment of g
                                                    {0., 0., 0.},
                                                    {0., 0., 0.}};
                         
                        for(int p=0; p < 3; ++p)
                        {
                            for(int q=0; q < 3; ++q)
                            {
                                eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + fluid1[i][j][k].p*velocity[p]*velocity[q];
                            }
                        }

                        double geq = fluid1[i][j][k].rho * fluid1[i][j][k].totalEnergy;
                        double gstr = fluid1[i][j][k].rho * fluid1[i][j][k].totalEnergy;

                        geq += dotproduct_Vec3(eq_heat_flux, velocity_set) / theta;
                        gstr += dotproduct_Vec3(str_heat_flux, velocity_set) / theta;

                        double matA[3][3] = {{eq_R_tensor[0][0]-fluid1[i][j][k].rho*fluid1[i][j][k].totalEnergy*theta, eq_R_tensor[0][1]                                                      , eq_R_tensor[0][2]},
                                             {eq_R_tensor[1][0]                                                      , eq_R_tensor[1][1]-fluid1[i][j][k].rho*fluid1[i][j][k].totalEnergy*theta, eq_R_tensor[1][2]},
                                             {eq_R_tensor[2][0]                                                      , eq_R_tensor[2][1]                                                      , eq_R_tensor[2][2]-fluid1[i][j][k].rho*fluid1[i][j][k].totalEnergy*theta}};
                        double matB[3][3] = {{cx[l]*cx[l]-theta  , cx[l]*cy[l]        , cx[l]*cz[l]},
                                             {cy[l]*cx[l]        , cy[l]*cy[l]-theta  , cy[l]*cz[l]},
                                             {cz[l]*cx[l]        , cz[l]*cy[l]        , cz[l]*cz[l]-theta}};
                        double result_AB = (matA[0][0]*matB[0][0] + matA[1][0]*matB[1][0] + matA[2][0]*matB[2][0]) + (matA[0][1]*matB[0][1] + matA[1][1]*matB[1][1] + matA[2][1]*matB[2][1]) + (matA[0][2]*matB[0][2] + matA[1][2]*matB[1][2] + matA[2][2]*matB[2][2]);
                        geq += result_AB/(2.0*theta*theta);
                        gstr += result_AB/(2.0*theta*theta);
                        
                        double weight = 1.0;
                        for (int m = 0; m < 3; ++m)
                        {
                            if (velocity_set[m] == 0) weight *= (1 - theta);
                            else weight *= theta / 2.0;
                        } 
                        geq = weight * geq;
                        gstr = weight * gstr;

                        double B;
                        double Z;
                        for (int m = 0; m < 3; ++m)
                        {
                            if (v_sqr(velocity_set[0], velocity_set[1], velocity_set[2]) == 0) B = 1;
                            else if (v_sqr(velocity_set[0], velocity_set[1], velocity_set[2]) == 1) B = -0.5*abs(velocity_set[m]);
                            else B = 0;
                            
                            Z = (1-3*theta)/(2*theta) * (eq_R_tensor[m][m]-theta*fluid1[i][j][k].rho*fluid1[i][j][k].totalEnergy);

                            geq += B*Z;
                            gstr += B*Z; 
                        }

                        fluid1[i][j][k].gpc[l] = fluid1[i][j][k].g[l] + fluid1[i][j][k].omega1*(geq-fluid1[i][j][k].g[l]) + (fluid1[i][j][k].omega-fluid1[i][j][k].omega1)*(gstr-fluid1[i][j][k].g[l]);
                    }
                    
                    #else
                        double uu = sq(fluid1[i][j][k].u) + sq(fluid1[i][j][k].v) + sq(fluid1[i][j][k].w);

                        for (int l = 0; l < npop; ++l)
                        {
                            double cu=cx[l]*fluid1[i][j][k].u+cy[l]*fluid1[i][j][k].v+cz[l]*fluid1[i][j][k].w;
                            double feq=wi[l]*fluid1[i][j][k].rho*(1.0+3.0*cu+4.5*cu*cu-1.5*uu);
                            fluid1[i][j][k].fpc[l]=(1.0-omega)*fluid1[i][j][k].f[l] + omega*feq;
                        }
                    #endif

                }
            }
        }
    } 
}

void LBM::Streaming()
{
    int i_nb, j_nb, k_nb;
    #pragma omp parallel for schedule(static, 1)
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k = 0; k<Nz; ++k)
            {
                if(fluid1[i][j][k].type==TYPE_F)
                {
                    for (int l=0; l < npop; ++l)
                    {
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        //---- Solid Boundary Condition ----------------------
                        if(fluid1[i_nb][j_nb][k_nb].type==TYPE_S)
                        {
                            fluid1[i][j][k].f[l] = fluid1[i][j][k].fpc[opposite[l]];
                            fluid1[i][j][k].g[l] = fluid1[i][j][k].gpc[opposite[l]];
                        }
                        //---- Inlet/Outlet Boundary Condition
                        else if (fluid1[i_nb][j_nb][k_nb].type==TYPE_E)
                        {
                            fluid1[i][j][k].f[l] = fluid1[i_nb][j_nb][k_nb].fpc[l];
                            fluid1[i][j][k].g[l] = fluid1[i_nb][j_nb][k_nb].gpc[l];
                        }
                        else //---- Periodic Boundary Condition --------------------
                        {
                            /*
                            if (i_nb < 1) i_nb = Nx-2;
                            else if(i_nb > Nx-2) i_nb = 1;

                            if (j_nb < 1) j_nb = Ny-2;
                            else if(j_nb > Ny-2) j_nb = 1;

                            if (k_nb < 1) k_nb = Nz-2;
                            else if(k_nb > Nz-2) k_nb = 1;

                            fluid1[i][j][k].f[l] = fluid1[i_nb][j_nb][k_nb].fpc[l];*/

                            
                            i_nb = ((i_nb - 1 + (Nx-2)) % (Nx-2)) + 1;
                            j_nb = ((j_nb - 1 + (Ny-2)) % (Ny-2)) + 1;
                            k_nb = ((k_nb - 1 + (Nz-2)) % (Nz-2)) + 1;
                            fluid1[i][j][k].f[l] = fluid1[i_nb][j_nb][k_nb].fpc[l];
                            fluid1[i][j][k].g[l] = fluid1[i_nb][j_nb][k_nb].gpc[l];
                        }
                    }
                }
            }
        }
    } 
}

void LBM::Quantity()
{
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int k = 0; k < Nz; ++k)
            {
                if (fluid1[i][j][k].type==TYPE_F)
                {
                    // Raw moment
                    double rho=0.0;
                    double rho_u=0.0;
                    double rho_v=0.0;
                    double rho_w=0.0;
                    
                    double rho_E=0.0;
                    double heat_flux_x=0.0;
                    double heat_flux_y=0.0;
                    double heat_flux_z=0.0;

                    for(int p=0; p < 3; ++p)
                    {
                        for(int q=0; q < 3; ++q)
                        {
                            fluid1[i][j][k].p_tensor[p][q] = 0.0;
                        }
                    }
                    
                    for (int l = 0; l < npop; ++l)
                    {
                        rho+=fluid1[i][j][k].f[l];
                        rho_u+=fluid1[i][j][k].f[l]*cx[l];
                        rho_v+=fluid1[i][j][k].f[l]*cy[l];
                        rho_w+=fluid1[i][j][k].f[l]*cz[l];
                        
                        double velocity_set[3] = {cx[l], cy[l], cz[l]};

                        for(int p=0; p < 3; ++p)
                        {
                            for(int q=0; q < 3; ++q)
                            {
                                fluid1[i][j][k].p_tensor[p][q] += fluid1[i][j][k].f[l]*velocity_set[p]*velocity_set[q];
                            }
                        }

                        rho_E+=fluid1[i][j][k].g[l];
                        heat_flux_x+=fluid1[i][j][k].g[l]*cx[l];
                        heat_flux_y+=fluid1[i][j][k].g[l]*cy[l];
                        heat_flux_z+=fluid1[i][j][k].g[l]*cz[l];
                    }

                    fluid1[i][j][k].rho=rho;
                    fluid1[i][j][k].u=rho_u/rho;
                    fluid1[i][j][k].v=rho_v/rho;
                    fluid1[i][j][k].w=rho_w/rho;

                    fluid1[i][j][k].totalEnergy=rho_E/rho;
                    fluid1[i][j][k].energy_flux[0]=heat_flux_x;
                    fluid1[i][j][k].energy_flux[1]=heat_flux_y;
                    fluid1[i][j][k].energy_flux[2]=heat_flux_z;
                    fluid1[i][j][k].internalEnergy=fluid1[i][j][k].totalEnergy - 0.5*v_sqr(fluid1[i][j][k].u,fluid1[i][j][k].v,fluid1[i][j][k].w);
                    double cv = fluid1[i][j][k].gas_const / (GAMMA - 1.0);
                    double cp = cv + fluid1[i][j][k].gas_const;
                    fluid1[i][j][k].temp = fluid1[i][j][k].internalEnergy / cv;
                    fluid1[i][j][k].enthalpy = cp*fluid1[i][j][k].temp;
                    fluid1[i][j][k].p = fluid1[i][j][k].rho*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp;

                    fluid1[i][j][k].mu = NU*fluid1[i][j][k].rho;
                    fluid1[i][j][k].omega = 2*fluid1[i][j][k].p*dt_sim / (fluid1[i][j][k].p*dt_sim + 2*fluid1[i][j][k].mu);
                    fluid1[i][j][k].conduc_coeff = fluid1[i][j][k].mu*cp/PR;
                    fluid1[i][j][k].omega1 = 2*fluid1[i][j][k].p*cp*dt_sim / (fluid1[i][j][k].p*cp*dt_sim + 2*fluid1[i][j][k].conduc_coeff);


                    // std::cout << "rho : " << fluid1[i][j][k].rho << std::endl; 
                    // std::cout << "R : " << fluid1[i][j][k].gas_const << std::endl;
                    // std::cout << "temp : " << fluid1[i][j][k].temp << std::endl;
                    
                    //std::cout << units.si_energy_mass(fluid1[i][j][k].internalEnergy) << std::endl;

                    // size_t rank = omp_get_thread_num();
                    // auto gas = sols[rank]->thermo();
                    // auto trans = sols[rank]->transport();

                    // gas->setState_UV(units.si_energy_mass(fluid1[i][j][k].internalEnergy), 1.0/units.si_rho(fluid1[i][j][k].rho));
                    // fluid1[i][j][k].mu =  units.mu(trans->viscosity());
                    // //fluid1[i][j][k].nu = NU; // debugg...
                    // fluid1[i][j][k].temp = units.temp(gas->temperature());
                    // //fluid1[i][j][k].temp = 2./3. * fluid1[i][j][k].internalEnergy / fluid1[i][j][k].gas_const; // debugg...
                    // fluid1[i][j][k].enthalpy = units.energy_mass(gas->enthalpy_mass());
                    // //fluid1[i][j][k].enthalpy = 5.0/2.0 * fluid1[i][j][k].gas_const * fluid1[i][j][k].temp; // debugg...
                    // fluid1[i][j][k].conduc_coeff = units.thermalConductivity(trans->thermalConductivity());
                    // fluid1[i][j][k].cp = units.cp(gas->cp_mass());
                    // fluid1[i][j][k].p=fluid1[i][j][k].rho*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp;

                    // fluid1[i][j][k].omega = 2*fluid1[i][j][k].p*dt_sim / (fluid1[i][j][k].p*dt_sim + 2*fluid1[i][j][k].mu);
                    // fluid1[i][j][k].omega1 = 2*fluid1[i][j][k].p*fluid1[i][j][k].cp*dt_sim / (fluid1[i][j][k].p*fluid1[i][j][k].cp*dt_sim + 2*fluid1[i][j][k].conduc_coeff);
                    // //double konst = PR * fluid1[i][j][k].omega / (2 - fluid1[i][j][k].omega); // debugg...
                    // //fluid1[i][j][k].omega1 = 2 * konst /(1 + konst); // debugg...
                    
                }
                else
                {
                    fluid1[i][j][k].rho=1.0;
                    fluid1[i][j][k].u=0.0;
                    fluid1[i][j][k].v=0.0;
                    fluid1[i][j][k].w=0.0;
                }
            }
        }
    }
    // std::cout << "omega  : " << fluid1[Nx/2][Ny/2][Nz/2].omega << std::endl;
    // std::cout << "omega1 : " << fluid1[Nx/2][Ny/2][Nz/2].omega1 << std::endl;
    // std::cout << "u : " << fluid1[Nx/2][Ny/2][Nz/2].u << std::endl;
    // std::cout << "p      : " << fluid1[Nx/2][Ny/2][Nz/2].p << std::endl;
    // std::cout << "lamda  : " << fluid1[Nx/2][Ny/2][Nz/2].conduc_coeff << std::endl;
}


// #include "./headers/lbm.hpp"
// #include <omp.h>
// #include "cantera.hpp"
// #include "lbm.hpp"
// #include "math_util.hpp"
// #include "units.hpp"
// #include <omp.h>
// #include <numeric>

// LBM::LBM(int Nx, int Ny, int Nz)
// {
//     this->Nx = Nx + 2;
//     this->Ny = Ny + 2;
//     this->Nz = Nz + 2;
//     // allocate memory for lattice
//     fluid1 = new LATTICE **[this->Nx];
//     for (int i = 0; i < this->Nx; ++i)
//     {
//         fluid1[i] = new LATTICE *[this->Ny];
//         for (int j = 0; j < this->Ny; ++j)
//         {
//             fluid1[i][j] = new LATTICE [this->Nz];
//         }
//     }
// }


// void LBM::Init()
// {
//     #pragma omp parallel for
//     for(int i = 0; i < Nx ; ++i)
//     {
//         for(int j = 0; j < Ny; ++j)
//         {
//             for(int k = 0; k < Nz; ++k)
//             {    
//                 if (fluid1[i][j][k].type == TYPE_F || fluid1[i][j][k].type == TYPE_E)     
//                 {    
//                     #ifdef LBM_EXTEND  
//                         double zeta = fluid1[i][j][k].gas_const * fluid1[i][j][k].temp;
//                         for (int l = 0; l < npop; ++l)
//                         {
//                             double feq = fluid1[i][j][k].rho;
//                             double P = 0.0;
                            
//                                 double eps = fluid1[i][j][k].u;
//                                 if (i > Nx - 3)
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].u,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].u,3)) - 4*(fluid1[i-1][j][k].rho*fluid1[i-1][j][k].u*(1-3*fluid1[i-1][j][k].gas_const*fluid1[i-1][j][k].temp)-fluid1[i-1][j][k].rho*pow(fluid1[i-1][j][k].u,3)) + (fluid1[i-2][j][k].rho*fluid1[i-2][j][k].u*(1-3*fluid1[i-2][j][k].gas_const*fluid1[i-2][j][k].temp)-fluid1[i-2][j][k].rho*pow(fluid1[i-2][j][k].u,3)))/(2*dx);
//                                 else
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].u,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(-3*(fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].u,3)) + 4*(fluid1[i+1][j][k].rho*fluid1[i+1][j][k].u*(1-3*fluid1[i+1][j][k].gas_const*fluid1[i+1][j][k].temp)-fluid1[i+1][j][k].rho*pow(fluid1[i+1][j][k].u,3)) - (fluid1[i+2][j][k].rho*fluid1[i+2][j][k].u*(1-3*fluid1[i+2][j][k].gas_const*fluid1[i+2][j][k].temp)-fluid1[i+2][j][k].rho*pow(fluid1[i+2][j][k].u,3)))/(2*dx);
                                
//                                 if (cx[l] == 0) feq *= (1 - P);
//                                 else if (cx[l] == 1) feq *= (eps+P)/2;
//                                 else if (cx[l] == -1) feq*= (-eps+P)/2;
                            
//                             #if NDIM == 2 || NDIM == 3
//                                 eps = fluid1[i][j][k].v;
//                                 if (j > Ny - 3)
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].v,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].v,3)) - 4*(fluid1[i][j-1][k].rho*fluid1[i][j-1][k].v*(1-3*fluid1[i][j-1][k].gas_const*fluid1[i][j-1][k].temp)-fluid1[i][j-1][k].rho*pow(fluid1[i][j-1][k].v,3)) + (fluid1[i][j-2][k].rho*fluid1[i][j-2][k].v*(1-3*fluid1[i][j-2][k].gas_const*fluid1[i][j-2][k].temp)-fluid1[i][j-2][k].rho*pow(fluid1[i][j-2][k].v,3)))/(2*dy);
//                                 else
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].v,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(-3*(fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].v,3)) + 4*(fluid1[i][j+1][k].rho*fluid1[i][j+1][k].v*(1-3*fluid1[i][j+1][k].gas_const*fluid1[i][j+1][k].temp)-fluid1[i][j+1][k].rho*pow(fluid1[i][j+1][k].v,3)) - (fluid1[i][j+2][k].rho*fluid1[i][j+2][k].v*(1-3*fluid1[i][j+2][k].gas_const*fluid1[i][j+2][k].temp)-fluid1[i][j+2][k].rho*pow(fluid1[i][j+2][k].v,3)))/(2*dy);
                                
//                                 if (cy[l] == 0) feq *= (1 - P);
//                                 else if (cy[l] == 1) feq *= (eps+P)/2;
//                                 else if (cy[l] == -1) feq*= (-eps+P)/2;
//                             #endif

//                             #if NDIM == 3
//                                 eps = fluid1[i][j][k].w;
//                                 if (k > Nz - 3)
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].w,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].w,3)) - 4*(fluid1[i][j][k-1].rho*fluid1[i][j][k-1].w*(1-3*fluid1[i][j][k-1].gas_const*fluid1[i][j][k-1].temp)-fluid1[i][j][k-1].rho*pow(fluid1[i][j][k-1].w,3)) + (fluid1[i][j][k-2].rho*fluid1[i][j][k-2].w*(1-3*fluid1[i][j][k-2].gas_const*fluid1[i][j][k-2].temp)-fluid1[i][j][k-2].rho*pow(fluid1[i][j][k-2].w,3)))/(2*dz);
//                                 else
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].w,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega) * (-3*(fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].w,3)) + 4*(fluid1[i][j][k+1].rho*fluid1[i][j][k+1].w*(1-3*fluid1[i][j][k+1].gas_const*fluid1[i][j][k+1].temp)-fluid1[i][j][k+1].rho*pow(fluid1[i][j][k+1].w,3)) - (fluid1[i][j][k+2].rho*fluid1[i][j][k+2].w*(1-3*fluid1[i][j][k+2].gas_const*fluid1[i][j][k+2].temp)-fluid1[i][j][k+2].rho*pow(fluid1[i][j][k+2].w,3)))/(2*dz);
                                
//                                 if (cz[l] == 0) feq *= (1 - P);
//                                 else if (cz[l] == 1) feq *= (eps+P)/2;
//                                 else if (cz[l] == -1) feq*= (-eps+P)/2;
//                             #endif
                            
//                             fluid1[i][j][k].fpc[l]=feq;
//                             fluid1[i][j][k].f[l]=feq;
//                         }

//                     #else
//                         double uu = fluid1[i][j][k].u*fluid1[i][j][k].u + fluid1[i][j][k].v*fluid1[i][j][k].v + fluid1[i][j][k].w*fluid1[i][j][k].w;

//                         for (int l = 0; l < npop; ++l)
//                         {
//                             double cu = cx[l]*fluid1[i][j][k].u + cy[l]*fluid1[i][j][k].v + cz[l]*fluid1[i][j][k].w;
//                             double feq = wi[l]*fluid1[i][j][k].rho*(1.0+3.0*cu+4.5*cu*cu-1.5*uu);
//                             fluid1[i][j][k].fpc[l]=feq;
//                             fluid1[i][j][k].f[l]=feq;
//                         }
//                     #endif
//                 }
//             }
//         }
//     }
// }

// void LBM::Collide_BGK()
// {
//     #pragma omp parallel for
//     for (int i = 0; i < Nx; ++i)
//     {
//         for (int j = 0; j < Ny; ++j)
//         {
//             for (int k = 0; k < Nz; ++k)
//             {
//                 if (fluid1[i][j][k].type==TYPE_F)
//                 {
                                       
//                     #ifdef LBM_EXTEND
//                         double zeta = fluid1[i][j][k].gas_const * fluid1[i][j][k].temp;
//                         for (int l = 0; l < npop; ++l)
//                         {
//                             double feq = fluid1[i][j][k].rho;
//                             double P = 0.0;
                            
//                                 double eps = fluid1[i][j][k].u;
//                                 if (i > Nx - 3)
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].u,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].u,3)) - 4*(fluid1[i-1][j][k].rho*fluid1[i-1][j][k].u*(1-3*fluid1[i-1][j][k].gas_const*fluid1[i-1][j][k].temp)-fluid1[i-1][j][k].rho*pow(fluid1[i-1][j][k].u,3)) + (fluid1[i-2][j][k].rho*fluid1[i-2][j][k].u*(1-3*fluid1[i-2][j][k].gas_const*fluid1[i-2][j][k].temp)-fluid1[i-2][j][k].rho*pow(fluid1[i-2][j][k].u,3)))/(2*dx);
//                                 else
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].u,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(-3*(fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].u,3)) + 4*(fluid1[i+1][j][k].rho*fluid1[i+1][j][k].u*(1-3*fluid1[i+1][j][k].gas_const*fluid1[i+1][j][k].temp)-fluid1[i+1][j][k].rho*pow(fluid1[i+1][j][k].u,3)) - (fluid1[i+2][j][k].rho*fluid1[i+2][j][k].u*(1-3*fluid1[i+2][j][k].gas_const*fluid1[i+2][j][k].temp)-fluid1[i+2][j][k].rho*pow(fluid1[i+2][j][k].u,3)))/(2*dx);
                                
//                                 if (cx[l] == 0) feq *= (1 - P);
//                                 else if (cx[l] == 1) feq *= (eps+P)/2;
//                                 else if (cx[l] == -1) feq*= (-eps+P)/2;
                            
//                             #if NDIM == 2 || NDIM == 3
//                                 eps = fluid1[i][j][k].v;
//                                 if (j > Ny - 3)
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].v,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].v,3)) - 4*(fluid1[i][j-1][k].rho*fluid1[i][j-1][k].v*(1-3*fluid1[i][j-1][k].gas_const*fluid1[i][j-1][k].temp)-fluid1[i][j-1][k].rho*pow(fluid1[i][j-1][k].v,3)) + (fluid1[i][j-2][k].rho*fluid1[i][j-2][k].v*(1-3*fluid1[i][j-2][k].gas_const*fluid1[i][j-2][k].temp)-fluid1[i][j-2][k].rho*pow(fluid1[i][j-2][k].v,3)))/(2*dy);
//                                 else
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].v,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(-3*(fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].v,3)) + 4*(fluid1[i][j+1][k].rho*fluid1[i][j+1][k].v*(1-3*fluid1[i][j+1][k].gas_const*fluid1[i][j+1][k].temp)-fluid1[i][j+1][k].rho*pow(fluid1[i][j+1][k].v,3)) - (fluid1[i][j+2][k].rho*fluid1[i][j+2][k].v*(1-3*fluid1[i][j+2][k].gas_const*fluid1[i][j+2][k].temp)-fluid1[i][j+2][k].rho*pow(fluid1[i][j+2][k].v,3)))/(2*dy);
                                
//                                 if (cy[l] == 0) feq *= (1 - P);
//                                 else if (cy[l] == 1) feq *= (eps+P)/2;
//                                 else if (cy[l] == -1) feq*= (-eps+P)/2;
//                             #endif

//                             #if NDIM == 3
//                                 eps = fluid1[i][j][k].w;
//                                 if (k > Nz - 3)
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].w,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].w,3)) - 4*(fluid1[i][j][k-1].rho*fluid1[i][j][k-1].w*(1-3*fluid1[i][j][k-1].gas_const*fluid1[i][j][k-1].temp)-fluid1[i][j][k-1].rho*pow(fluid1[i][j][k-1].w,3)) + (fluid1[i][j][k-2].rho*fluid1[i][j][k-2].w*(1-3*fluid1[i][j][k-2].gas_const*fluid1[i][j][k-2].temp)-fluid1[i][j][k-2].rho*pow(fluid1[i][j][k-2].w,3)))/(2*dz);
//                                 else
//                                     P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].w,2) + dt_sim*(2-fluid1[i][j][k].omega)/(2*fluid1[i][j][k].rho*fluid1[i][j][k].omega) * (-3*(fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].w,3)) + 4*(fluid1[i][j][k+1].rho*fluid1[i][j][k+1].w*(1-3*fluid1[i][j][k+1].gas_const*fluid1[i][j][k+1].temp)-fluid1[i][j][k+1].rho*pow(fluid1[i][j][k+1].w,3)) - (fluid1[i][j][k+2].rho*fluid1[i][j][k+2].w*(1-3*fluid1[i][j][k+2].gas_const*fluid1[i][j][k+2].temp)-fluid1[i][j][k+2].rho*pow(fluid1[i][j][k+2].w,3)))/(2*dz);
                                
//                                 if (cz[l] == 0) feq *= (1 - P);
//                                 else if (cz[l] == 1) feq *= (eps+P)/2;
//                                 else if (cz[l] == -1) feq*= (-eps+P)/2;
//                             #endif
                            
//                             fluid1[i][j][k].fpc[l]=(1.0-fluid1[i][j][k].omega)*fluid1[i][j][k].f[l]+fluid1[i][j][k].omega*feq;
//                         }

//                     #else
//                         double uu = u*u + v*v + w*w;

//                         for (int l = 0; l < npop; ++l)
//                         {
//                             double cu=cx[l]*u+cy[l]*v+cz[l]*w;
//                             double feq=wi[l]*rho*(1.0+3.0*cu+4.5*cu*cu-1.5*uu);
//                             fluid1[i][j][k].fpc[l]=(1.0-fluid1[i][j][k].omega)*fluid1[i][j][k].f[l] + fluid1[i][j][k].omega*feq;
//                         }
//                     #endif

//                 }
//             }
//         }
//     } 
// }

// void LBM::Streaming()
// {
//     int i_nb, j_nb, k_nb;
//     #pragma omp parallel for
//     for(int i=0; i<Nx; ++i)
//     {
//         for(int j=0; j<Ny; ++j)
//         {
//             for(int k = 0; k<Nz; ++k)
//             {
//                 if(fluid1[i][j][k].type==TYPE_F)
//                 {
//                     for (int l=0; l < npop; ++l)
//                     {
//                         i_nb = i - cx[l];
//                         j_nb = j - cy[l];
//                         k_nb = k - cz[l];

//                         //---- Solid Boundary Condition ----------------------
//                         if(fluid1[i_nb][j_nb][k_nb].type==TYPE_S)
//                         {
//                             fluid1[i][j][k].f[l] = fluid1[i][j][k].fpc[opposite[l]];
//                         }
//                         //---- Inlet/Outlet Boundary Condition
//                         else if (fluid1[i_nb][j_nb][k_nb].type==TYPE_E)
//                         {
//                             fluid1[i][j][k].f[l] = fluid1[i_nb][j_nb][k_nb].fpc[l];
//                         }
//                         else //---- Periodic Boundary Condition --------------------
//                         {
//                             /*
//                             if (i_nb < 1) i_nb = Nx-2;
//                             else if(i_nb > Nx-2) i_nb = 1;

//                             if (j_nb < 1) j_nb = Ny-2;
//                             else if(j_nb > Ny-2) j_nb = 1;

//                             if (k_nb < 1) k_nb = Nz-2;
//                             else if(k_nb > Nz-2) k_nb = 1;

//                             fluid1[i][j][k].f[l] = fluid1[i_nb][j_nb][k_nb].fpc[l];*/

                            
//                             i_nb = ((i_nb - 1 + (Nx-2)) % (Nx-2)) + 1;
//                             j_nb = ((j_nb - 1 + (Ny-2)) % (Ny-2)) + 1;
//                             k_nb = ((k_nb - 1 + (Nz-2)) % (Nz-2)) + 1;
//                             fluid1[i][j][k].f[l] = fluid1[i_nb][j_nb][k_nb].fpc[l];
//                         }
//                     }
//                 }
//             }
//         }
//     } 
// }

// void LBM::Quantity()
// {
//     #pragma omp parallel for
//     for (int i = 0; i < Nx; ++i)
//     {
//         for (int j = 0; j < Ny; ++j)
//         {
//             for (int k = 0; k < Nz; ++k)
//             {
//                 if (fluid1[i][j][k].type==TYPE_F)
//                 {
//                     // Raw moment
//                     double rho=0.0;
//                     double rho_u=0.0;
//                     double rho_v=0.0;
//                     double rho_w=0.0;
//                     for (int l = 0; l < npop; ++l)
//                     {
//                         rho+=fluid1[i][j][k].f[l];
//                         rho_u+=fluid1[i][j][k].f[l]*cx[l];
//                         rho_v+=fluid1[i][j][k].f[l]*cy[l];
//                         rho_w+=fluid1[i][j][k].f[l]*cz[l];
//                     }
//                     fluid1[i][j][k].rho=rho;
//                     fluid1[i][j][k].u=rho_u/rho;
//                     fluid1[i][j][k].v=rho_v/rho;
//                     fluid1[i][j][k].w=rho_w/rho;
//                 }
//                 else
//                 {
//                     fluid1[i][j][k].rho=1.0;
//                     fluid1[i][j][k].u=0.0;
//                     fluid1[i][j][k].v=0.0;
//                     fluid1[i][j][k].w=0.0;
//                 }
//             }
//         }
//     }
// }
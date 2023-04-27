
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
    mixture = new LATTICE **[this->Nx];
    for (int i = 0; i < this->Nx; ++i)
    {
        mixture[i] = new LATTICE *[this->Ny];
        for (int j = 0; j < this->Ny; ++j)
        {
            mixture[i][j] = new LATTICE [this->Nz];
        }
    }

    // int nThreads  = omp_get_max_threads();

    // for(int i = 0; i < nThreads; ++i)
    // {
    //     auto sol = newSolution("gri30.yaml", "gri30");
    //     sols.emplace_back(sol);
    // }
}

LBM::LBM(int Nx, int Ny, int Nz, int nSpecies)
{
    this->Nx = Nx + 2;
    this->Ny = Ny + 2;
    this->Nz = Nz + 2;
    this->nSpecies = nSpecies;

    // allocate memory for mixture
    mixture = new LATTICE **[this->Nx];
    for (int i = 0; i < this->Nx; ++i)
    {
        mixture[i] = new LATTICE *[this->Ny];
        for (int j = 0; j < this->Ny; ++j)
        {
            mixture[i][j] = new LATTICE [this->Nz];
        }
    }

    // allocate memory for species
    this->species.resize(nSpecies);
    for(int p = 0; p < nSpecies; ++p)
    {
        species[p] = new SPECIES **[this->Nx];
        for (int i = 0; i < this->Nx; ++i)
        {
            species[p][i] = new SPECIES *[this->Ny];
            for (int j = 0; j < this->Ny; ++j)
            {
                species[p][i][j] = new SPECIES [this->Nz];
            }
        }
    }
}

void LBM::calculate_moment(LATTICE ***fluid)
{
    std::vector<std::vector<std::vector<double>>> dQdevx(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevy(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevz(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));

    // #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int k = 0; k < Nz; ++k)
            {
                if (fluid[i][j][k].type==TYPE_F)
                {
                    // Raw moment
                    double rho = 0.0;
                    double rhou = 0.0;
                    double rhov = 0.0;
                    double rhow = 0.0;
                    double rhoe = 0.0;

                    double heat_flux_x=0.0;
                    double heat_flux_y=0.0;
                    double heat_flux_z=0.0;

                    for(int p=0; p < 3; ++p)
                    {
                        for(int q=0; q < 3; ++q)
                        {
                            fluid[i][j][k].p_tensor[p][q] = 0.0;
                        }
                    }
                    
                    for (int l = 0; l < npop; ++l)
                    {
                        rho+=fluid[i][j][k].f[l];
                        rhou+=fluid[i][j][k].f[l]*cx[l];
                        rhov+=fluid[i][j][k].f[l]*cy[l];
                        rhow+=fluid[i][j][k].f[l]*cz[l];
                        
                        double velocity_set[3] = {cx[l], cy[l], cz[l]};

                        for(int p=0; p < 3; ++p)
                        {
                            for(int q=0; q < 3; ++q)
                            {
                                fluid[i][j][k].p_tensor[p][q] += fluid[i][j][k].f[l]*velocity_set[p]*velocity_set[q];
                            }
                        }

                        rhoe += fluid[i][j][k].g[l];
                        heat_flux_x += fluid[i][j][k].g[l]*cx[l];
                        heat_flux_y += fluid[i][j][k].g[l]*cy[l];
                        heat_flux_z += fluid[i][j][k].g[l]*cz[l];
                    }

                    fluid[i][j][k].rho = rho;
                    fluid[i][j][k].rhou = rhou;
                    fluid[i][j][k].rhov = rhov;
                    fluid[i][j][k].rhow = rhow;
                    fluid[i][j][k].rhoe = rhoe;
                    fluid[i][j][k].energy_flux[0]=heat_flux_x;
                    fluid[i][j][k].energy_flux[1]=heat_flux_y;
                    fluid[i][j][k].energy_flux[2]=heat_flux_z;
                    
                    // std::cout << "rho : " << fluid[i][j][k].rho << std::endl; 
                    // std::cout << "R : " << fluid[i][j][k].gas_const << std::endl;
                    // std::cout << "temp : " << fluid[i][j][k].temp << std::endl;
                    
                    //std::cout << units.si_energy_mass(fluid[i][j][k].internalEnergy) << std::endl;

                    // size_t rank = omp_get_thread_num();
                    // auto gas = sols[rank]->thermo();
                    // auto trans = sols[rank]->transport();

                    // gas->setState_UV(units.si_energy_mass(fluid[i][j][k].internalEnergy), 1.0/units.si_rho(fluid[i][j][k].rho));
                    // fluid[i][j][k].mu =  units.mu(trans->viscosity());
                    // //fluid[i][j][k].nu = NU; // debugg...
                    // fluid[i][j][k].temp = units.temp(gas->temperature());
                    // //fluid[i][j][k].temp = 2./3. * fluid[i][j][k].internalEnergy / fluid[i][j][k].gas_const; // debugg...
                    // fluid[i][j][k].enthalpy = units.energy_mass(gas->enthalpy_mass());
                    // //fluid[i][j][k].enthalpy = 5.0/2.0 * fluid[i][j][k].gas_const * fluid[i][j][k].temp; // debugg...
                    // fluid[i][j][k].conduc_coeff = units.thermalConductivity(trans->thermalConductivity());
                    // fluid[i][j][k].cp = units.cp(gas->cp_mass());
                    // fluid[i][j][k].p=fluid[i][j][k].rho*fluid[i][j][k].gas_const*fluid[i][j][k].temp;

                    // fluid[i][j][k].omega = 2*fluid[i][j][k].p*dt_sim / (fluid[i][j][k].p*dt_sim + 2*fluid[i][j][k].mu);
                    // fluid[i][j][k].omega1 = 2*fluid[i][j][k].p*fluid[i][j][k].cp*dt_sim / (fluid[i][j][k].p*fluid[i][j][k].cp*dt_sim + 2*fluid[i][j][k].conduc_coeff);
                    // //double konst = PR * fluid[i][j][k].omega / (2 - fluid[i][j][k].omega); // debugg...
                    // //fluid[i][j][k].omega1 = 2 * konst /(1 + konst); // debugg...

                    double velocity[3] = {  fluid[i][j][k].rhou/fluid[i][j][k].rho,
                                            fluid[i][j][k].rhov/fluid[i][j][k].rho, 
                                            fluid[i][j][k].rhow/fluid[i][j][k].rho};
                    double internalEnergy=fluid[i][j][k].rhoe/fluid[i][j][k].rho - 0.5*v_sqr(velocity[0], velocity[1], velocity[2]);
                    double gas_const = GAS_CONST/fluid[i][j][k].mmass;
                    double cv = gas_const / (GAMMA - 1.0);
                    fluid[i][j][k].temp = internalEnergy / cv;
                    fluid[i][j][k].p = fluid[i][j][k].rho*gas_const*fluid[i][j][k].temp;
                }
                else
                {
                    fluid[i][j][k].rho=1.0;
                    fluid[i][j][k].rhou=0.0;
                    fluid[i][j][k].rhov=0.0;
                    fluid[i][j][k].rhow=0.0;
                    fluid[i][j][k].rhoe=0.0;
                }
            }
        }
    }

    #pragma region calculate Qdev
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (fluid[i][j][k].type == TYPE_F || fluid[i][j][k].type == TYPE_E)     
                {
                    double gas_const = GAS_CONST/fluid[i][j][k].mmass;
                    if (i == 1)
                    {
                        dQdevx[i][j][k] = ((fluid[i+1][j][k].rhou*(1-3*gas_const*fluid[i+1][j][k].temp)-fluid[i+1][j][k].rho*cb(fluid[i+1][j][k].rhou/fluid[i+1][j][k].rho)) - (fluid[i][j][k].rhou*(1-3*gas_const*fluid[i][j][k].temp)-fluid[i][j][k].rho*cb(fluid[i][j][k].rhou/fluid[i][j][k].rho)) ) /(dx);
                    }
                    else
                    {
                        dQdevx[i][j][k] = ((fluid[i][j][k].rhou*(1-3*gas_const*fluid[i][j][k].temp)-fluid[i][j][k].rho*cb(fluid[i][j][k].rhou/fluid[i][j][k].rho)) - (fluid[i-1][j][k].rhou*(1-3*gas_const*fluid[i-1][j][k].temp)-fluid[i-1][j][k].rho*cb(fluid[i-1][j][k].rhou/fluid[i-1][j][k].rho)) ) /(dx);
                    }
                    if (j == 1)
                    {
                        dQdevy[i][j][k] = ((fluid[i][j+1][k].rhov*(1-3*gas_const*fluid[i][j+1][k].temp)-fluid[i][j+1][k].rho*cb(fluid[i][j+1][k].rhov/fluid[i][j+1][k].rho)) - (fluid[i][j][k].rhov*(1-3*gas_const*fluid[i][j][k].temp)-fluid[i][j][k].rho*cb(fluid[i][j][k].rhov/fluid[i][j][k].rho))) /(dy);
                    }
                    else
                    {
                        dQdevy[i][j][k] = ((fluid[i][j][k].rhov*(1-3*gas_const*fluid[i][j][k].temp)-fluid[i][j][k].rho*cb(fluid[i][j][k].rhov/fluid[i][j][k].rho)) - (fluid[i][j-1][k].rhov*(1-3*gas_const*fluid[i][j-1][k].temp)-fluid[i][j-1][k].rho*cb(fluid[i][j-1][k].rhov/fluid[i][j-1][k].rho)) ) /(dy);
                    }
                    if (k == 1)
                    {
                        dQdevz[i][j][k] = ((fluid[i][j][k+1].rhow*(1-3*gas_const*fluid[i][j][k+1].temp)-fluid[i][j][k+1].rho*cb(fluid[i][j][k+1].rhow/fluid[i][j][k+1].rho)) - (fluid[i][j][k].rhow*(1-3*gas_const*fluid[i][j][k].temp)-fluid[i][j][k].rho*cb(fluid[i][j][k].rhow/fluid[i][j][k].rho)))  /(dz);
                    }
                    else
                    {
                        dQdevz[i][j][k] = ((fluid[i][j][k].rhow*(1-3*gas_const*fluid[i][j][k].temp)-fluid[i][j][k].rho*cb(fluid[i][j][k].rhow/fluid[i][j][k].rho)) - (fluid[i][j][k-1].rhow*(1-3*gas_const*fluid[i][j][k-1].temp)-fluid[i][j][k-1].rho*cb(fluid[i][j][k-1].rhow/fluid[i][j][k-1].rho)) ) /(dz);
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
                if (fluid[i][j][k].type == TYPE_F || fluid[i][j][k].type == TYPE_E)     
                {
                    if (i == Nx-2) fluid[i][j][k].dQdevx = limiterVanleer(dQdevx[i-1][j][k],dQdevx[i][j][k]);
                    else fluid[i][j][k].dQdevx = limiterVanleer(dQdevx[i][j][k],dQdevx[i+1][j][k]);
                    
                    //if (fluid[i][j][k].dQdevx != 0.0) std::cout << "dQdev : " << fluid[i][j][k].dQdevx << " | " << i << " | " << j << " | " << k << std::endl;

                    if (j == Ny-2) fluid[i][j][k].dQdevy = limiterVanleer(dQdevy[i][j-1][k],dQdevy[i][j][k]);
                    else fluid[i][j][k].dQdevy = limiterVanleer(dQdevy[i][j][k],dQdevy[i][j+1][k]);

                    if (k == Nz-2) fluid[i][j][k].dQdevz = limiterVanleer(dQdevz[i][j][k-1],dQdevz[i][j][k]);
                    else fluid[i][j][k].dQdevz = limiterVanleer(dQdevz[i][j][k],dQdevz[i][j][k+1]);
                }
            }
        }
    }
    #pragma endregion

}

double LBM::calculate_feq(int l, double rho, double velocity[], double theta, double corr[])
{
    double P = 0.0;
    double eps = 0.0;
    double feq = rho;

        eps = velocity[0];
        P = theta + sq(velocity[0]) + corr[0];
        if (cx[l] == 0) feq *= (1 - P);
        else if (cx[l] == 1) feq *= (eps+P)/2;
        else if (cx[l] == -1) feq*= (-eps+P)/2;
    
    #if NDIM == 2 || NDIM == 3
        eps = velocity[1];
        P = theta + sq(velocity[1]) + corr[1];;
        if (cy[l] == 0) feq *= (1 - P);
        else if (cy[l] == 1) feq *= (eps+P)/2;
        else if (cy[l] == -1) feq*= (-eps+P)/2;
    #endif

    #if NDIM == 3
        eps = velocity[2];
        P = theta + sq(velocity[2]) + corr[2];
        if (cz[l] == 0) feq *= (1 - P);
        else if (cz[l] == 1) feq *= (eps+P)/2;
        else if (cz[l] == -1) feq*= (-eps+P)/2;
    #endif

    return feq;
}


double LBM::calculate_geq(int l, double rhoe, double eq_heat_flux[], double eq_R_tensor[][3], double theta)
{
    double geq = rhoe;
    
    double velocity_set[3] = {cx[l], cy[l], cz[l]};
    geq += dotproduct_Vec3(eq_heat_flux, velocity_set) / theta;

    double matA[3][3] = {{eq_R_tensor[0][0]-rhoe*theta, eq_R_tensor[0][1]           , eq_R_tensor[0][2]},
                         {eq_R_tensor[1][0]           , eq_R_tensor[1][1]-rhoe*theta, eq_R_tensor[1][2]},
                         {eq_R_tensor[2][0]           , eq_R_tensor[2][1]           , eq_R_tensor[2][2]-rhoe*theta}};
    double matB[3][3] = {{cx[l]*cx[l]-theta  , cx[l]*cy[l]        , cx[l]*cz[l]},
                         {cy[l]*cx[l]        , cy[l]*cy[l]-theta  , cy[l]*cz[l]},
                         {cz[l]*cx[l]        , cz[l]*cy[l]        , cz[l]*cz[l]-theta}};
    double result_AB = (matA[0][0]*matB[0][0] + matA[1][0]*matB[1][0] + matA[2][0]*matB[2][0]) + (matA[0][1]*matB[0][1] + matA[1][1]*matB[1][1] + matA[2][1]*matB[2][1]) + (matA[0][2]*matB[0][2] + matA[1][2]*matB[1][2] + matA[2][2]*matB[2][2]);
    geq += result_AB/(2.0*theta*theta);
    
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
        
        Z = (1-3*theta)/(2*theta) * (eq_R_tensor[m][m]-theta*rhoe);

        geq += B*Z;
    }

    return geq;
}


void LBM::Init(LATTICE ***fluid)
{
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (fluid[i][j][k].type == TYPE_F || fluid[i][j][k].type == TYPE_E)     
                {                
                    double velocity[3] = {  fluid[i][j][k].rhou/fluid[i][j][k].rho,
                                            fluid[i][j][k].rhov/fluid[i][j][k].rho, 
                                            fluid[i][j][k].rhow/fluid[i][j][k].rho};
                    double gas_const = GAS_CONST/fluid[i][j][k].mmass;
                    double cv = gas_const / (GAMMA - 1.0);
                    double cp = cv + gas_const;
                    double internal_energy = cv * fluid[i][j][k].temp;
                    fluid[i][j][k].p = fluid[i][j][k].rho*gas_const*fluid[i][j][k].temp;
                    fluid[i][j][k].rhoe = fluid[i][j][k].rho*(internal_energy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]));
                    double theta = gas_const*fluid[i][j][k].temp;   
                    double enthalpy = cp * fluid[i][j][k].temp; // H = Cp * T = (Cv + 1) * T
                    double total_enthalpy = enthalpy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    double eq_heat_flux[3] = {  total_enthalpy*fluid[i][j][k].rhou,
                                                total_enthalpy*fluid[i][j][k].rhov,
                                                total_enthalpy*fluid[i][j][k].rhow};
                    double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                {0., 0., 0.},
                                                {0., 0., 0.}};
                    double eq_R_tensor[3][3] = {{0., 0., 0.},    // second-order moment of g
                                                {0., 0., 0.},
                                                {0., 0., 0.}};

                    for(int p=0; p < 3; ++p)
                    {
                            for(int q=0; q < 3; ++q)
                            {
                                eq_p_tensor[p][q] = (p==q) ? fluid[i][j][k].p+fluid[i][j][k].rho*velocity[p]*velocity[q] : fluid[i][j][k].rho*velocity[p]*velocity[q]; 
                                eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + fluid[i][j][k].p*velocity[p]*velocity[q];
                            }  
                    }

                    double corr[3] = {0, 0, 0}; 

                    for (int l = 0; l < npop; ++l)
                    {
                        // ------------- Mass and Momentum Initialization -----------------------------
                        fluid[i][j][k].f[l]=calculate_feq(l, fluid[i][j][k].rho, velocity, theta, corr);
                        fluid[i][j][k].g[l]=calculate_geq(l, fluid[i][j][k].rhoe, eq_heat_flux, eq_R_tensor, theta);
                    }     
                }
            }
        }
    }
}

void LBM::Collide(LATTICE ***fluid)
{
    calculate_moment(fluid);
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (fluid[i][j][k].type == TYPE_F || fluid[i][j][k].type == TYPE_E)     
                {   
                    double velocity[3] = {  fluid[i][j][k].rhou/fluid[i][j][k].rho,
                                            fluid[i][j][k].rhov/fluid[i][j][k].rho, 
                                            fluid[i][j][k].rhow/fluid[i][j][k].rho};
                    double gas_const = GAS_CONST/fluid[i][j][k].mmass;
                    double cv = gas_const / (GAMMA - 1.0);
                    double cp = cv + gas_const;
                    double theta = gas_const*fluid[i][j][k].temp;
                    double mu = NU*fluid[i][j][k].rho;
                    double omega = 2*fluid[i][j][k].p*dt_sim / (fluid[i][j][k].p*dt_sim + 2*mu);
                    double conduc_coeff = mu*cp/PR;
                    double omega1 = 2*fluid[i][j][k].p*cp*dt_sim / (fluid[i][j][k].p*cp*dt_sim + 2*conduc_coeff);
                    double enthalpy = cp * fluid[i][j][k].temp; // H = Cp * T = (Cv + 1) * T
                    double total_enthalpy = enthalpy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    double eq_heat_flux[3] = {  total_enthalpy*fluid[i][j][k].rhou,
                                                total_enthalpy*fluid[i][j][k].rhov,
                                                total_enthalpy*fluid[i][j][k].rhow};
                    double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                {0., 0., 0.},
                                                {0., 0., 0.}};
                    double eq_R_tensor[3][3] = {{0., 0., 0.},    // second-order moment of g
                                                {0., 0., 0.},
                                                {0., 0., 0.}};

                    for(int p=0; p < 3; ++p)
                    {
                            for(int q=0; q < 3; ++q)
                            {
                                eq_p_tensor[p][q] = (p==q) ? fluid[i][j][k].p+fluid[i][j][k].rho*velocity[p]*velocity[q] : fluid[i][j][k].rho*velocity[p]*velocity[q]; 
                                eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + fluid[i][j][k].p*velocity[p]*velocity[q];
                            }  
                    }

                    double str_heat_flux[3] ={  fluid[i][j][k].energy_flux[0] - velocity[0]*(fluid[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) - velocity[1]*(fluid[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) - velocity[2]*(fluid[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) - 0.5*dt_sim*velocity[0]*fluid[i][j][k].dQdevx,
                                                fluid[i][j][k].energy_flux[1] - velocity[0]*(fluid[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) - velocity[1]*(fluid[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) - velocity[2]*(fluid[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) - 0.5*dt_sim*velocity[1]*fluid[i][j][k].dQdevy,
                                                fluid[i][j][k].energy_flux[2] - velocity[0]*(fluid[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) - velocity[1]*(fluid[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) - velocity[2]*(fluid[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) - 0.5*dt_sim*velocity[2]*fluid[i][j][k].dQdevz};

                    double corr[3] = {  dt_sim*(2-omega)/(2*fluid[i][j][k].rho*omega)*fluid[i][j][k].dQdevx,
                                        dt_sim*(2-omega)/(2*fluid[i][j][k].rho*omega)*fluid[i][j][k].dQdevy,
                                        dt_sim*(2-omega)/(2*fluid[i][j][k].rho*omega)*fluid[i][j][k].dQdevz};

                    // double corr[3] = {0, 0, 0};

                    for (int l = 0; l < npop; ++l)
                    {
                        // ------------- Mass and Momentum Initialization -----------------------------
                        double feq = calculate_feq(l, fluid[i][j][k].rho, velocity, theta, corr);
                        double geq = calculate_geq(l, fluid[i][j][k].rhoe, eq_heat_flux, eq_R_tensor, theta);
                        double gstr = calculate_geq(l, fluid[i][j][k].rhoe, str_heat_flux, eq_R_tensor, theta);

                        fluid[i][j][k].fpc[l]= (1.0-omega)*fluid[i][j][k].f[l] + omega*feq;
                        fluid[i][j][k].gpc[l] = fluid[i][j][k].g[l] + omega1*(geq-fluid[i][j][k].g[l]) + (omega-omega1)*(gstr-fluid[i][j][k].g[l]);
                    }     
                }
            }
        }
    }
}


void LBM::Streaming(LATTICE ***fluid)
{
    int i_nb, j_nb, k_nb;
    #pragma omp parallel for schedule(static, 1)
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k = 0; k<Nz; ++k)
            {
                if(fluid[i][j][k].type==TYPE_F)
                {
                    for (int l=0; l < npop; ++l)
                    {
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        //---- Solid Boundary Condition ----------------------
                        if(fluid[i_nb][j_nb][k_nb].type==TYPE_S)
                        {
                            fluid[i][j][k].f[l] = fluid[i][j][k].fpc[opposite[l]];
                            fluid[i][j][k].g[l] = fluid[i][j][k].gpc[opposite[l]];
                        }
                        //---- Inlet/Outlet Boundary Condition
                        else if (fluid[i_nb][j_nb][k_nb].type==TYPE_E)
                        {
                            fluid[i][j][k].f[l] = fluid[i_nb][j_nb][k_nb].fpc[l];
                            fluid[i][j][k].g[l] = fluid[i_nb][j_nb][k_nb].gpc[l];
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

                            fluid[i][j][k].f[l] = fluid[i_nb][j_nb][k_nb].fpc[l];*/

                            
                            i_nb = ((i_nb - 1 + (Nx-2)) % (Nx-2)) + 1;
                            j_nb = ((j_nb - 1 + (Ny-2)) % (Ny-2)) + 1;
                            k_nb = ((k_nb - 1 + (Nz-2)) % (Nz-2)) + 1;
                            fluid[i][j][k].f[l] = fluid[i_nb][j_nb][k_nb].fpc[l];
                            fluid[i][j][k].g[l] = fluid[i_nb][j_nb][k_nb].gpc[l];
                        }
                    }
                }
            }
        }
    } 
}


// void LBM::Init(SPECIES*** species)
// {
//     #pragma omp parallel for schedule(static, 1)
//     for(int i = 0; i < Nx ; ++i)
//     {
//         for(int j = 0; j < Ny; ++j)
//         {
//             for(int k = 0; k < Nz; ++k)
//             {    
//                 if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E)     
//                 {                
//                     double velocity[3] = {  species[i][j][k].rhou/species[i][j][k].rho,
//                                             species[i][j][k].rhov/species[i][j][k].rho, 
//                                             species[i][j][k].rhow/species[i][j][k].rho};
//                     double gas_const = GAS_CONST/mixture[i][j][k].mmass;
//                     double cv = gas_const / (GAMMA - 1.0);
//                     double cp = cv + gas_const;
//                     double internal_energy = cv * fluid[i][j][k].temp;
//                     fluid[i][j][k].p = fluid[i][j][k].rho*gas_const*fluid[i][j][k].temp;
//                     fluid[i][j][k].rhoe = fluid[i][j][k].rho*(internal_energy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]));
//                     double theta = gas_const*fluid[i][j][k].temp;   

//                     double corr[3] = {0, 0, 0}; 

//                     for (int l = 0; l < npop; ++l)
//                     {
//                         // ------------- Mass and Momentum Initialization -----------------------------
//                         fluid[i][j][k].f[l]=calculate_feq(l, fluid[i][j][k].rho, velocity, theta, corr);
//                     }     
//                 }
//             }
//         }
//     }
// }




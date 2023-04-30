
#include "cantera.hpp"
#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include <omp.h>
#include <numeric>
#include <vector>


LBM::LBM(int Nx, int Ny, int Nz, double NU)
{
    this->Nx = Nx + 2;
    this->Ny = Ny + 2;
    this->Nz = Nz + 2;
    this->nu = NU;

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
}

#ifdef MULTICOMP
LBM::LBM(int Nx, int Ny, int Nz, std::vector<std::string> species)
{
    this->Nx = Nx + 2;
    this->Ny = Ny + 2;
    this->Nz = Nz + 2;
    this->speciesName = species;
    this->nSpecies = species.size();

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
        this->species[p] = new SPECIES **[this->Nx];
        for (int i = 0; i < this->Nx; ++i)
        {
            this->species[p][i] = new SPECIES *[this->Ny];
            for (int j = 0; j < this->Ny; ++j)
            {
                this->species[p][i][j] = new SPECIES [this->Nz];
            }
        }
    }
}
#endif

void LBM::calculate_moment()
{
    std::vector<std::vector<std::vector<double>>> dQdevx(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevy(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevz(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));

    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int k = 0; k < Nz; ++k)
            {
                if (mixture[i][j][k].type==TYPE_F)
                {
                    // mixture moment
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
                            mixture[i][j][k].p_tensor[p][q] = 0.0;
                        }
                    }
                    
                    for (int l = 0; l < npop; ++l)
                    {
                        rho+=mixture[i][j][k].f[l];
                        rhou+=mixture[i][j][k].f[l]*cx[l];
                        rhov+=mixture[i][j][k].f[l]*cy[l];
                        rhow+=mixture[i][j][k].f[l]*cz[l];
                        
                        double velocity_set[3] = {cx[l], cy[l], cz[l]};

                        for(int p=0; p < 3; ++p)
                        {
                            for(int q=0; q < 3; ++q)
                            {
                                mixture[i][j][k].p_tensor[p][q] += mixture[i][j][k].f[l]*velocity_set[p]*velocity_set[q];
                            }
                        }

                        rhoe += mixture[i][j][k].g[l];
                        heat_flux_x += mixture[i][j][k].g[l]*cx[l];
                        heat_flux_y += mixture[i][j][k].g[l]*cy[l];
                        heat_flux_z += mixture[i][j][k].g[l]*cz[l];
                    }

                    mixture[i][j][k].rho = rho;
                    mixture[i][j][k].u = rhou / rho;
                    mixture[i][j][k].v = rhov / rho;
                    mixture[i][j][k].w = rhow / rho;
                    mixture[i][j][k].rhoe = rhoe;
                    mixture[i][j][k].energy_flux[0]=heat_flux_x;
                    mixture[i][j][k].energy_flux[1]=heat_flux_y;
                    mixture[i][j][k].energy_flux[2]=heat_flux_z;

                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};
                    double internalEnergy=mixture[i][j][k].rhoe/mixture[i][j][k].rho - 0.5*v_sqr(velocity[0], velocity[1], velocity[2]);
                    #ifndef MULTICOMP
                        double cv = gas_const / (gamma - 1.0);
                        mixture[i][j][k].temp = internalEnergy / cv;
                        mixture[i][j][k].p = mixture[i][j][k].rho*gas_const*mixture[i][j][k].temp;
                    #else
                        // species moment
                        auto gas = mixture[i][j][k].sol->thermo();
                        std::vector<double> Y (gas->nSpecies());
                        std::vector<double> X (gas->nSpecies());

                        for(int a = 0; a < nSpecies; ++a)
                        {
                            double rho = 0.0;
                            double rhou = 0.0;
                            double rhov = 0.0;
                            double rhow = 0.0;
                            
                            for (int l = 0; l < npop; ++l)
                            {
                                rho  += species[a][i][j][k].f[l];
                                rhou += species[a][i][j][k].f[l]*cx[l];
                                rhov += species[a][i][j][k].f[l]*cy[l];
                                rhow += species[a][i][j][k].f[l]*cz[l];
                            }

                            species[a][i][j][k].rho = rho;
                            species[a][i][j][k].u = rhou / rho;
                            species[a][i][j][k].v = rhov / rho;
                            species[a][i][j][k].w = rhow / rho;
                            Y[gas->speciesIndex(speciesName[a])]  = species[a][i][j][k].rho / mixture[i][j][k].rho;
                        } 
                        gas->setState_UV(units.si_energy_mass(internalEnergy), 1/units.si_rho(rho));
                        gas->setMassFractions(&Y[0]);
                        gas->getMoleFractions(&X[0]);
                        for(int a = 0; a < nSpecies; ++a) species[a][i][j][k].X = X[gas->speciesIndex(speciesName[a])];
                        mixture[i][j][k].temp = units.temp(gas->temperature());
                        mixture[i][j][k].p = units.p(gas->pressure());
                    #endif

                }
                else
                {
                    mixture[i][j][k].rho=1.0;
                    mixture[i][j][k].u=0.0;
                    mixture[i][j][k].v=0.0;
                    mixture[i][j][k].w=0.0;
                    mixture[i][j][k].rhoe=0.0;
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
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E)     
                {
                    #ifdef MULTICOMP
                        auto gas = mixture[i][j][k].sol->thermo();
                        gas_const = units.cp(Cantera::GasConstant/gas->meanMolecularWeight());
                    #endif
                    if (i == 1)
                    {
                        dQdevx[i][j][k] = ( (mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1-3*gas_const*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u)) - (mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u)) ) /(dx);
                    }
                    else
                    {
                        dQdevx[i][j][k] = ( (mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u)) - (mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1-3*gas_const*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u)) ) /(dx);
                    }
                    if (j == 1)
                    {
                        dQdevy[i][j][k] = ( (mixture[i][j+1][k].rho*mixture[i][j+1][k].v*(1-3*gas_const*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].v)) - (mixture[i][j][k].rho*mixture[i][j][k].v*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v)) ) /(dy);
                    }
                    else
                    {
                        dQdevy[i][j][k] = ( (mixture[i][j][k].rho*mixture[i][j][k].v*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v)) - (mixture[i][j-1][k].rho*mixture[i][j-1][k].v*(1-3*gas_const*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].v)) ) /(dy);
                    }
                    if (k == 1)
                    {
                        dQdevz[i][j][k] = ( (mixture[i][j][k+1].rho*mixture[i][j][k+1].w*(1-3*gas_const*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].w)) - (mixture[i][j][k].rho*mixture[i][j][k].w*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w)) ) /(dz);
                    }
                    else
                    {
                        dQdevz[i][j][k] = ( (mixture[i][j][k].rho*mixture[i][j][k].w*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w))- (mixture[i][j][k-1].rho*mixture[i][j][k-1].w*(1-3*gas_const*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].w)) ) /(dz);
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
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E)     
                {
                    if (i == Nx-2) mixture[i][j][k].dQdevx = limiterVanleer(dQdevx[i-1][j][k],dQdevx[i][j][k]);
                    else mixture[i][j][k].dQdevx = limiterVanleer(dQdevx[i][j][k],dQdevx[i+1][j][k]);
                    
                    //if (mixture[i][j][k].dQdevx != 0.0) std::cout << "dQdev : " << mixture[i][j][k].dQdevx << " | " << i << " | " << j << " | " << k << std::endl;

                    if (j == Ny-2) mixture[i][j][k].dQdevy = limiterVanleer(dQdevy[i][j-1][k],dQdevy[i][j][k]);
                    else mixture[i][j][k].dQdevy = limiterVanleer(dQdevy[i][j][k],dQdevy[i][j+1][k]);

                    if (k == Nz-2) mixture[i][j][k].dQdevz = limiterVanleer(dQdevz[i][j][k-1],dQdevz[i][j][k]);
                    else mixture[i][j][k].dQdevz = limiterVanleer(dQdevz[i][j][k],dQdevz[i][j][k+1]);
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


void LBM::Init()
{
    //#pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                std::cout << i << " " << j << " " << k << std::endl;
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E)     
                {            
                    #ifndef MULTICOMP    
                        double velocity[3] = {  mixture[i][j][k].u,
                                                mixture[i][j][k].v, 
                                                mixture[i][j][k].w};
                        double cv = gas_const / (gamma - 1.0);
                        double cp = cv + gas_const;
                        double internal_energy = cv * mixture[i][j][k].temp;
                        mixture[i][j][k].p = mixture[i][j][k].rho*gas_const*mixture[i][j][k].temp;
                        mixture[i][j][k].rhoe = mixture[i][j][k].rho*(internal_energy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]));
                        double theta = gas_const*mixture[i][j][k].temp;   
                        double enthalpy = cp * mixture[i][j][k].temp; // H = Cp * T = (Cv + 1) * T
                    #else
                        mixture[i][j][k].sol = Cantera::newSolution("gri30.yaml", "gri30");
                        auto gas = mixture[i][j][k].sol->thermo();
                        std::vector <double> X (gas->nSpecies());
                        for(int a = 0; a < nSpecies; ++a) X[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].X;
                        gas->setState_TPX(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].p), &X[0]);
                       
                        double internal_energy = units.energy_mass(gas->intEnergy_mass());
                        double enthalpy = units.energy_mass(gas->enthalpy_mass());
                        double velocity[3] = {  mixture[i][j][k].u,
                                                mixture[i][j][k].v, 
                                                mixture[i][j][k].w};
                        mixture[i][j][k].rho = units.rho(gas->density());
                        mixture[i][j][k].rhoe = mixture[i][j][k].rho*(internal_energy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]));
                        double theta = units.energy_mass(gas->RT()/gas->meanMolecularWeight());

                        X.clear();
                    #endif
                        

                    double total_enthalpy = enthalpy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    double eq_heat_flux[3] = {  total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].u,
                                                total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].v,
                                                total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].w};
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
                                eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 
                                eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + mixture[i][j][k].p*velocity[p]*velocity[q];
                            }  
                    }

                    double corr[3] = {0, 0, 0}; 

                    for (int l = 0; l < npop; ++l)
                    {
                        // ------------- Mass and Momentum Initialization -----------------------------
                        mixture[i][j][k].f[l]=calculate_feq(l, mixture[i][j][k].rho, velocity, theta, corr);
                        mixture[i][j][k].g[l]=calculate_geq(l, mixture[i][j][k].rhoe, eq_heat_flux, eq_R_tensor, theta);
                    }     

                    #ifdef MULTICOMP
                    for(int a = 0; a < nSpecies; ++a)
                    {
                        double velocity[3] = {  species[a][i][j][k].u,
                                                species[a][i][j][k].v, 
                                                species[a][i][j][k].w};
                        double theta = units.energy_mass(gas->RT() / gas->molecularWeight(gas->speciesIndex(speciesName[a])));

                        for (int l = 0; l < npop; ++l)
                        {
                            species[a][i][j][k].f[l]=calculate_feq(l, species[a][i][j][k].rho, velocity, theta, corr);
                        }   
                    }
                    #endif
                }
            }
        }
    }
}

void LBM::Collide()
{
    calculate_moment();
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E)     
                {   
                    #ifndef MULTICOMP
                        double cv = gas_const / (gamma - 1.0);
                        double cp = cv + gas_const;
                        double theta = gas_const*mixture[i][j][k].temp;
                        double mu = nu*mixture[i][j][k].rho;
                        double conduc_coeff = mu*cp/prtl;
                        double enthalpy = cp * mixture[i][j][k].temp; // H = Cp * T = (Cv + 1) * T
                    
                    #else
                        auto gas =  mixture[i][j][k].sol->thermo();
                        double cp = units.cp(gas->cp_mass());
                        double enthalpy = units.energy_mass(gas->enthalpy_mass());

                        auto trans = mixture[i][j][k].sol->transport();
                        double mu = units.mu(trans->viscosity());
                        double conduc_coeff = units.thermalConductivity(trans->thermalConductivity());

                        double theta = units.energy_mass(gas->RT()/gas->meanMolecularWeight());
                    #endif

                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};
                    double omega1 = 2*mixture[i][j][k].p*cp*dt_sim / (mixture[i][j][k].p*cp*dt_sim + 2*conduc_coeff);
                    double omega = 2*mixture[i][j][k].p*dt_sim / (mixture[i][j][k].p*dt_sim + 2*mu);
                    double total_enthalpy = enthalpy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    double eq_heat_flux[3] = {  total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].u,
                                                total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].v,
                                                total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].w};
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
                                eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 
                                eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + mixture[i][j][k].p*velocity[p]*velocity[q];
                            }  
                    }

                    //double q_diff [3] += {(omega1/(omega-omega1)) * mixture[i][j][k].rho * };

                    double str_heat_flux[3] ={  mixture[i][j][k].energy_flux[0] - velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) - velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) - velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) - 0.5*dt_sim*velocity[0]*mixture[i][j][k].dQdevx,
                                                mixture[i][j][k].energy_flux[1] - velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) - velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) - velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) - 0.5*dt_sim*velocity[1]*mixture[i][j][k].dQdevy,
                                                mixture[i][j][k].energy_flux[2] - velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) - velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) - velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) - 0.5*dt_sim*velocity[2]*mixture[i][j][k].dQdevz};

                    double corr[3] = {  dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*mixture[i][j][k].dQdevx,
                                        dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*mixture[i][j][k].dQdevy,
                                        dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*mixture[i][j][k].dQdevz};

                    for (int l = 0; l < npop; ++l)
                    {
                        // ------------- Mass and Momentum Initialization -----------------------------
                        double feq = calculate_feq(l, mixture[i][j][k].rho, velocity, theta, corr);
                        double geq = calculate_geq(l, mixture[i][j][k].rhoe, eq_heat_flux, eq_R_tensor, theta);
                        double gstr = calculate_geq(l, mixture[i][j][k].rhoe, str_heat_flux, eq_R_tensor, theta);

                        mixture[i][j][k].fpc[l]= (1.0-omega)*mixture[i][j][k].f[l] + omega*feq;
                        mixture[i][j][k].gpc[l] = mixture[i][j][k].g[l] + omega1*(geq-mixture[i][j][k].g[l]) + (omega-omega1)*(gstr-mixture[i][j][k].g[l]);
                    }     
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
                if(mixture[i][j][k].type==TYPE_F)
                {
                    for (int l=0; l < npop; ++l)
                    {
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        //---- Solid Boundary Condition ----------------------
                        if(mixture[i_nb][j_nb][k_nb].type==TYPE_S)
                        {
                            mixture[i][j][k].f[l] = mixture[i][j][k].fpc[opposite[l]];
                            mixture[i][j][k].g[l] = mixture[i][j][k].gpc[opposite[l]];
                        }
                        //---- Inlet/Outlet Boundary Condition
                        else if (mixture[i_nb][j_nb][k_nb].type==TYPE_E)
                        {
                            mixture[i][j][k].f[l] = mixture[i_nb][j_nb][k_nb].fpc[l];
                            mixture[i][j][k].g[l] = mixture[i_nb][j_nb][k_nb].gpc[l];
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

                            mixture[i][j][k].f[l] = mixture[i_nb][j_nb][k_nb].fpc[l];*/

                            
                            i_nb = ((i_nb - 1 + (Nx-2)) % (Nx-2)) + 1;
                            j_nb = ((j_nb - 1 + (Ny-2)) % (Ny-2)) + 1;
                            k_nb = ((k_nb - 1 + (Nz-2)) % (Nz-2)) + 1;
                            mixture[i][j][k].f[l] = mixture[i_nb][j_nb][k_nb].fpc[l];
                            mixture[i][j][k].g[l] = mixture[i_nb][j_nb][k_nb].gpc[l];
                        }
                    }
                }
            }
        }
    } 
}

#ifdef MULTICOMP
void LBM::Collide_Species()
{
    calculate_moment();
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E)     
                {   
                    double mass_frac[nSpecies];         
                    double tau_ab[nSpecies][nSpecies];
                    double D_ab[nSpecies][nSpecies];
                    double invtau_a[nSpecies];
                    double *du[nSpecies];
                    double *dv[nSpecies];
                    double *dw[nSpecies];
                    double mmass[nSpecies]; 
                    Eigen::MatrixXd mat_A(nSpecies, nSpecies);
                    Eigen::VectorXd vec_b(nSpecies);
                    Eigen::VectorXd vec_u(nSpecies);
                    Eigen::VectorXd vec_v(nSpecies);
                    Eigen::VectorXd vec_w(nSpecies);
                    Eigen::VectorXd sol(nSpecies);

                    auto gas = mixture[i][j][k].sol->thermo();                   
                    auto trans = mixture[i][j][k].sol->transport();

                    for(int a = 0; a < nSpecies; ++a) 
                    {
                        mass_frac[a] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                        mmass[a] = units.M(gas->molecularWeight(gas->speciesIndex(speciesName[a])));
                    }
                    double mix_mmass = units.M(gas->meanMolecularWeight());
                    double theta = units.energy(gas->RT());
                    int ld = gas->nSpecies();
                    double d[ld * ld];
                    trans->getBinaryDiffCoeffs(ld, d);

                    for(int a = 0; a < nSpecies; ++a)
                    {                       
                        for(int b = a; b < nSpecies; ++b)
                        {
                            D_ab[a][b] = d[ld*b + a];
                            tau_ab[a][b] = (mmass[a]*mmass[b]/(mix_mmass*theta))*D_ab[a][b];
                        }
                        
                        for(int b = 0; b < nSpecies; ++b)
                            if (a != b)
                            {
                                if (a < b) invtau_a[a] += mass_frac[b] / tau_ab[a][b]; 
                                else invtau_a[a] += mass_frac[b] / tau_ab[b][a];   
                            }

                        for(int b = 0; b < nSpecies; ++b)
                            if (a == b) mat_A(a, b) = 1 + dt_sim * invtau_a[a] / 2.0;
                            else 
                                if (a < b) mat_A(a, b) = - dt_sim / 2 * mass_frac[b] / tau_ab[a][b]; 
                                else mat_A(a, b) = - dt_sim / 2 * mass_frac[b] / tau_ab[b][a]; 

                        vec_u(a) = species[a][i][j][k].u - mixture[i][j][k].u;
                        vec_v(a) = species[a][i][j][k].v - mixture[i][j][k].v;
                        vec_w(a) = species[a][i][j][k].w - mixture[i][j][k].w;
                    }  
                        
                    sol = mat_A.colPivHouseholderQr().solve(vec_u);
                    *du = sol.data();
                    sol = mat_A.colPivHouseholderQr().solve(vec_v);
                    *dv = sol.data();
                    sol = mat_A.colPivHouseholderQr().solve(vec_w);
                    *dw = sol.data();                  

                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};

                    for (int l = 0; l < npop; ++l)
                    {
                        double feq[nSpecies];
                        double fstr[nSpecies];

                        for (int a = 0; a < nSpecies; ++a)
                        {
                            double velocity_spec[3] = { velocity[0] + *du[a],
                                                        velocity[1] + *dv[a], 
                                                        velocity[2] + *dw[a]};
                            species[a][i][j][k].u = velocity_spec[0];
                            species[a][i][j][k].v = velocity_spec[1];
                            species[a][i][j][k].w = velocity_spec[2];
                            auto gas = mixture[i][j][k].sol->thermo();
                            double theta = units.energy_mass( gas->RT()/gas->molecularWeight(gas->speciesIndex(speciesName[a])) );
                            double corr[3] = {0, 0, 0};
                            feq[a] = calculate_feq(l, species[a][i][j][k].rho, velocity, theta, corr);
                            fstr[a] = calculate_feq(l, species[a][i][j][k].rho, velocity_spec, theta, corr);
                        }

                        for(int a = 0; a < nSpecies; ++a)
                        {
                            double F_a = 0.0;
                            double beta = 0.0;
                            for(int b = 0; b < nSpecies; ++b)
                                if (a !=b)
                                {
                                    if (a < b)
                                        F_a += mass_frac[a]/tau_ab[a][b]*(feq[b]-fstr[b]);
                                    else 
                                        F_a += mass_frac[a]/tau_ab[b][a]*(feq[b]-fstr[b]);
                                }
                            
                            beta = dt_sim / (2*(1/invtau_a[a]) + dt_sim);
                            
                            species[a][i][j][k].fpc[l] = species[a][i][j][k].f[l] + 2.0*beta*(feq[a]-species[a][i][j][k].f[l]) + dt_sim*(beta-1)*F_a;
                        }
                    }  

                }
            }
        }
    }
}
#endif

void LBM::run(int nstep, int tout)
{
    this->nstep = nstep;
    this->tout = tout;

    std::cout << "-- Setup Done --" << std::endl;

    // initialize the distribution function 
    Init();  
    std::cout << "-- Initialization Done --" << std::endl;

    calculate_moment();

    // initialize time step & Save the macroscopic at t=0
    int step = 0;
    OutputVTK(step, this);
    OutputKeEns(step, this);
    
    // Simulation loop
    for (step = 1; step <= nstep; ++step)
    {
        Collide();   // collision step
        // std::cout << "-- Collision Done --" << std::endl;
        Streaming();     // streaming step & BC
        // std::cout << "-- Streaming Done --" << std::endl;

        if (step % tout == 0)
        {
            //std::cout << "Step : " << step << std::endl;
            OutputKeEns(step, this);
            //OutputVTK(step, lb); // Save the macroscopic quantity
        }

    }
}



#include "cantera.hpp"
#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include <omp.h>
#include <numeric>
#include <vector>

#ifdef MULTICOMP
    std::vector<std::shared_ptr<Cantera::Solution>> sols;
#endif

LBM::LBM(int Nx, int Ny, int Nz, double nu){
    // Number lattice used in simulation
    this->Nx = Nx + 2;  // + 2 for ghost lattice in the boundary
    this->Ny = Ny + 2;
    this->Nz = Nz + 2;
    this->nu = nu;

    // allocate memory for lattice
    mixture = new MIXTURE **[this->Nx];
    for (int i = 0; i < this->Nx; ++i){
        mixture[i] = new MIXTURE *[this->Ny];
        for (int j = 0; j < this->Ny; ++j)
            mixture[i][j] = new MIXTURE [this->Nz];
    }
}

#ifdef MULTICOMP
LBM::LBM(int Nx, int Ny, int Nz, std::vector<std::string> species){
    // Number lattice used in simulation
    this->Nx = Nx + 2;  // + 2 for ghost lattice in the boundary
    this->Ny = Ny + 2;
    this->Nz = Nz + 2;
    this->speciesName = species;
    this->nSpecies = species.size();

    // allocate memory for mixture
    mixture = new MIXTURE **[this->Nx];
    for (int i = 0; i < this->Nx; ++i){
        mixture[i] = new MIXTURE *[this->Ny];
        for (int j = 0; j < this->Ny; ++j){
            mixture[i][j] = new MIXTURE [this->Nz];
        }
    }

    // allocate memory for species
    this->species.resize(nSpecies);
    for(size_t p = 0; p < nSpecies; ++p){
        this->species[p] = new SPECIES **[this->Nx];
        for (int i = 0; i < this->Nx; ++i){
            this->species[p][i] = new SPECIES *[this->Ny];
            for (int j = 0; j < this->Ny; ++j){
                this->species[p][i][j] = new SPECIES [this->Nz];
            }
        }
    }

    // Create Cantera's Solution object
    int nThreads = omp_get_max_threads();
    std::cout << "nThreads : " << nThreads << std::endl;
    for(int i = 0; i < nThreads; ++i){
        // auto sol = Cantera::newSolution("gri30.yaml", "gri30", "multicomponent");
        auto sol = Cantera::newSolution("h2o2.yaml", "ohmech","mixture-averaged");
        sols.push_back(sol);
    }

}
#endif

void LBM::calculate_moment()
{
    std::vector<std::vector<std::vector<double>>> dQdevx(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevy(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevz(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<std::vector<double>>>> delYx(nSpecies, std::vector<std::vector<std::vector<double>>> (Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz))));
    std::vector<std::vector<std::vector<std::vector<double>>>> delYy(nSpecies, std::vector<std::vector<std::vector<double>>> (Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz))));
    std::vector<std::vector<std::vector<std::vector<double>>>> delYz(nSpecies, std::vector<std::vector<std::vector<double>>> (Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz))));


    #ifdef PARALLEL 
    #pragma omp parallel for schedule(static, 1) 
    #endif
    for (int i = 0; i < Nx; ++i){
        #ifdef MULTICOMP
        int rank = omp_get_thread_num();
        #endif
        for (int j = 0; j < Ny; ++j){
            for (int k = 0; k < Nz; ++k){
                if (mixture[i][j][k].type==TYPE_F){
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
                        for(int q=0; q < 3; ++q)
                            mixture[i][j][k].p_tensor[p][q] = 0.0;
                    
                    for (int l = 0; l < npop; ++l){
                        rho+=mixture[i][j][k].f[l];
                        rhou+=mixture[i][j][k].f[l]*cx[l];
                        rhov+=mixture[i][j][k].f[l]*cy[l];
                        rhow+=mixture[i][j][k].f[l]*cz[l];
                        
                        double velocity_set[3] = {cx[l], cy[l], cz[l]};
                        for(int p=0; p < 3; ++p)
                            for(int q=0; q < 3; ++q)
                                mixture[i][j][k].p_tensor[p][q] += mixture[i][j][k].f[l]*velocity_set[p]*velocity_set[q];

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
                    mixture[i][j][k].energy_flux[0] = heat_flux_x;
                    mixture[i][j][k].energy_flux[1] = heat_flux_y;
                    mixture[i][j][k].energy_flux[2] = heat_flux_z;
                                        
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
                        auto gas = sols[rank]->thermo();
                        std::vector<double> Y (gas->nSpecies());
                        std::vector<double> X (gas->nSpecies());
                        
                        for(size_t a = 0; a < nSpecies; ++a){

                        #ifndef FD
                            double rho = 0.0;
                            double rhou = 0.0;
                            double rhov = 0.0;
                            double rhow = 0.0;
                            
                            for (int l = 0; l < npop; ++l){
                                rho  += species[a][i][j][k].f[l];
                                rhou += species[a][i][j][k].f[l]*cx[l];
                                rhov += species[a][i][j][k].f[l]*cy[l];
                                rhow += species[a][i][j][k].f[l]*cz[l];
                            }

                            species[a][i][j][k].rho = rho;

                            if(rho != 0){
                                species[a][i][j][k].u = rhou / rho;
                                species[a][i][j][k].v = rhov / rho;
                                species[a][i][j][k].w = rhow / rho;
                            }
                            else{
                                species[a][i][j][k].u = 0.0;
                                species[a][i][j][k].v = 0.0;
                                species[a][i][j][k].w = 0.0;
                            }
                            Y[gas->speciesIndex(speciesName[a])]  = species[a][i][j][k].rho / mixture[i][j][k].rho;
                            species[a][i][j][k].rho_dot = 0.0;

                        #elif defined FD
                            species[a][i][j][k].rho = species[a][i][j][k].rho_n; 
                            Y[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho; 
                        #endif
                        }

                        gas->setMassFractions(&Y[0]);
                        gas->setState_UV(units.si_energy_mass(internalEnergy), 1.0/units.si_rho(mixture[i][j][k].rho),1.0e-12);
                        gas->getMoleFractions(&X[0]);
                        for(size_t a = 0; a < nSpecies; ++a) species[a][i][j][k].X = X[gas->speciesIndex(speciesName[a])];

                        mixture[i][j][k].temp = units.temp(gas->temperature()); 
                        mixture[i][j][k].p = units.p(gas->pressure()); 
                        gas_const = units.cp(Cantera::GasConstant/gas->meanMolecularWeight());

                    #endif
                }
            }
        }
    }

    fill_BC();

    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (mixture[i][j][k].type == TYPE_F|| mixture[i][j][k].type == TYPE_E)     
                {
                    #ifdef MULTICOMP
                        int rank = omp_get_thread_num();
                        auto gas = sols[rank]->thermo();
                        std::vector<double> Y (gas->nSpecies());
                        for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho ;
                        gas->setMassFractions(&Y[0]);
                        gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));
                        
                        gas_const = units.cp(Cantera::GasConstant/gas->meanMolecularWeight());
                    #endif

                    if (i == 1)
                    {
                        dQdevx[i][j][k] = ( (mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1-3*gas_const*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u)) - (mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u)) ) /(dx);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            delYx[a][i][j][k] = ((species[a][i+1][j][k].rho / mixture[i+1][j][k].rho) - (species[a][i][j][k].rho / mixture[i][j][k].rho)) / (dx);
                        #endif
                    }
                    else
                    {
                        dQdevx[i][j][k] = ( (mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u)) - (mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1-3*gas_const*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u)) ) /(dx);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            delYx[a][i][j][k] = ((species[a][i][j][k].rho / mixture[i][j][k].rho) - (species[a][i-1][j][k].rho / mixture[i-1][j][k].rho)) / (dx);
                        #endif
                    }

                    if (j == 1)
                    {
                        dQdevy[i][j][k] = ( (mixture[i][j+1][k].rho*mixture[i][j+1][k].v*(1-3*gas_const*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].v)) - (mixture[i][j][k].rho*mixture[i][j][k].v*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v)) ) /(dy);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            delYy[a][i][j][k] = ((species[a][i][j+1][k].rho / mixture[i][j+1][k].rho) - (species[a][i][j][k].rho / mixture[i][j][k].rho)) / (dy);
                        #endif
                    }
                    else
                    {
                        dQdevy[i][j][k] = ( (mixture[i][j][k].rho*mixture[i][j][k].v*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v)) - (mixture[i][j-1][k].rho*mixture[i][j-1][k].v*(1-3*gas_const*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].v)) ) /(dy);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            delYy[a][i][j][k] = ((species[a][i][j][k].rho / mixture[i][j][k].rho) - (species[a][i][j-1][k].rho / mixture[i][j-1][k].rho)) / (dy);
                        #endif
                    }

                    if (k == 1)
                    {
                        dQdevz[i][j][k] = ( (mixture[i][j][k+1].rho*mixture[i][j][k+1].w*(1-3*gas_const*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].w)) - (mixture[i][j][k].rho*mixture[i][j][k].w*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w)) ) /(dz);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            delYz[a][i][j][k] = ((species[a][i][j][k+1].rho / mixture[i][j][k+1].rho) - (species[a][i][j][k].rho / mixture[i][j][k].rho)) / (dz);
                        #endif
                    }
                    else
                    {
                        dQdevz[i][j][k] = ( (mixture[i][j][k].rho*mixture[i][j][k].w*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w))- (mixture[i][j][k-1].rho*mixture[i][j][k-1].w*(1-3*gas_const*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].w)) ) /(dz);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            delYz[a][i][j][k] = ((species[a][i][j][k].rho / mixture[i][j][k].rho) - (species[a][i][j][k-1].rho / mixture[i][j][k-1].rho)) / (dz);
                        #endif
                    }                
                    
                    #ifdef MULTICOMP
                        // for(size_t a = 0; a < nSpecies; a++)
                        // {                                                     
                        //     species[a][i][j][k].delYx = ((species[a][i+1][j][k].rho / mixture[i+1][j][k].rho) - (species[a][i-1][j][k].rho / mixture[i-1][j][k].rho)) / (2*dx);
                        //     species[a][i][j][k].delYy = ((species[a][i][j+1][k].rho / mixture[i][j+1][k].rho) - (species[a][i][j-1][k].rho / mixture[i][j-1][k].rho)) / (2*dy);
                        //     species[a][i][j][k].delYz = ((species[a][i][j][k+1].rho / mixture[i][j][k+1].rho) - (species[a][i][j][k-1].rho / mixture[i][j][k-1].rho)) / (2*dz);
                        //     // std::cout << i << " | " << idx_p << " | " << idx_n << " | " << species[a][i][j][k].delYx << " | " << species[a][i][j][k].delYy << " | " << species[a][i][j][k].delYz << std::endl;
                        // }
                    #endif

                }
            }
        }
    }

    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (mixture[i][j][k].type == TYPE_F|| mixture[i][j][k].type == TYPE_E)     
                {
                    if (i == Nx-2) 
                    {
                        mixture[i][j][k].dQdevx = LIMITER_TYPE(dQdevx[i-1][j][k], dQdevx[i][j][k]);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k].delYx = LIMITER_TYPE(delYx[a][i-1][j][k], delYx[a][i][j][k]);
                        #endif
                    }
                    else 
                    {
                        mixture[i][j][k].dQdevx = LIMITER_TYPE(dQdevx[i][j][k],dQdevx[i+1][j][k]);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k].delYx = LIMITER_TYPE(delYx[a][i][j][k], delYx[a][i+1][j][k]);
                        #endif  
                    }

                    if (j == Ny-2) 
                    {
                        mixture[i][j][k].dQdevy = LIMITER_TYPE(dQdevy[i][j-1][k],dQdevy[i][j][k]);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k].delYy = LIMITER_TYPE(delYy[a][i][j-1][k], delYy[a][i][j][k]);
                        #endif
                    }
                    else 
                    {
                        mixture[i][j][k].dQdevy = LIMITER_TYPE(dQdevy[i][j][k],dQdevy[i][j+1][k]);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k].delYy = LIMITER_TYPE(delYy[a][i][j][k], delYy[a][i][j+1][k]);
                        #endif
                    }

                    if (k == Nz-2) 
                    {
                        mixture[i][j][k].dQdevz = LIMITER_TYPE(dQdevz[i][j][k-1],dQdevz[i][j][k]);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k].delYz = LIMITER_TYPE(delYz[a][i][j][k-1], delYz[a][i][j][k]);
                        #endif
                    }
                    else 
                    {
                        mixture[i][j][k].dQdevz = LIMITER_TYPE(dQdevz[i][j][k],dQdevz[i][j][k+1]);
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k].delYz = LIMITER_TYPE(delYz[a][i][j][k], delYz[a][i][j][k+1]);
                        #endif
                    }
                }
            }
        }
    }
}

double LBM::calculate_feq(int l, double rho, double velocity[], double theta, double corr[]){
    double eps = 0.0;
    double P = 0.0;
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


double LBM::calculate_geq(int l, double rhoe, double eq_heat_flux[], double eq_R_tensor[][3], double theta){
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
    for (int m = 0; m < 3; ++m){
        if (velocity_set[m] == 0) weight *= (1 - theta);
        else weight *= theta / 2.0;
    } 
    geq = weight * geq;

    double B;
    double Z;
    for (int m = 0; m < 3; ++m){
        if (v_sqr(velocity_set[0], velocity_set[1], velocity_set[2]) == 0) B = 1;
        else if (v_sqr(velocity_set[0], velocity_set[1], velocity_set[2]) == 1) B = -0.5*abs(velocity_set[m]);
        else B = 0;
        
        Z = (1-3*theta)/(2*theta) * (eq_R_tensor[m][m]-theta*rhoe);

        geq += B*Z;
    }

    return geq;
}

// Initialize the disribution function from known macroscopic properties
void LBM::Init(){
    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){    
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E){            
                    #ifndef MULTICOMP    
                        double velocity[3] = {  mixture[i][j][k].u,
                                                mixture[i][j][k].v, 
                                                mixture[i][j][k].w};
                        double cv = gas_const / (gamma - 1.0);
                        double cp = cv + gas_const;
                        double internal_energy = cv * mixture[i][j][k].temp;
                        mixture[i][j][k].rho = mixture[i][j][k].p / (gas_const*mixture[i][j][k].temp);
                        mixture[i][j][k].rhoe = mixture[i][j][k].rho*(internal_energy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]));
                        double theta = gas_const*mixture[i][j][k].temp;   
                        double enthalpy = cp * mixture[i][j][k].temp; // H = Cp * T = (Cv + 1) * T
                    #else
                        // Initiate Cantera object
                        int rank = omp_get_thread_num();
                        auto gas = sols[rank]->thermo();
                        std::vector <double> X (gas->nSpecies());
                        for(size_t a = 0; a < nSpecies; ++a) {
                            X[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].X;
                        }
                        gas->setMoleFractions(&X[0]);
                        gas->setState_TP(units.si_temp(mixture[i][j][k].temp), units.si_p(mixture[i][j][k].p));

                        // Calculate other macropscopic properties [pressure, internal energy, enthalpy, total energy]
                        mixture[i][j][k].rho = units.rho(gas->density());
                        double internal_energy = units.energy_mass(gas->intEnergy_mass());
                        double enthalpy = units.energy_mass(gas->enthalpy_mass());
                        double velocity[3] = {  mixture[i][j][k].u,
                                                mixture[i][j][k].v, 
                                                mixture[i][j][k].w};
                        mixture[i][j][k].rhoe = mixture[i][j][k].rho*(internal_energy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]));
                        double theta = units.energy_mass(gas->RT()/gas->meanMolecularWeight());
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

                    for(int p=0; p < 3; ++p){
                        for(int q=0; q < 3; ++q){
                            eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 
                            eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + mixture[i][j][k].p*velocity[p]*velocity[q];
                        }  
                    }

                    double corr[3] = {0, 0, 0}; 

                    for (int l = 0; l < npop; ++l){
                        // ------------- Mass and Momentum Distribution Function Initialization -----------------------------
                        mixture[i][j][k].f[l]=calculate_feq(l, mixture[i][j][k].rho, velocity, theta, corr);
                        mixture[i][j][k].g[l]=calculate_geq(l, mixture[i][j][k].rhoe, eq_heat_flux, eq_R_tensor, theta);
                    }     

                    #if defined MULTICOMP && !defined FD
                        // Species distribution function Initialization
                        for(size_t a = 0; a < nSpecies; ++a){
                            species[a][i][j][k].rho = gas->massFraction(gas->speciesIndex(speciesName[a])) * mixture[i][j][k].rho;
                            theta = units.energy_mass(gas->RT() / gas->molecularWeight(gas->speciesIndex(speciesName[a])) );

                            for (int l = 0; l < npop; ++l)
                                species[a][i][j][k].f[l]=calculate_feq(l, species[a][i][j][k].rho, velocity, theta, corr);
                        }

                    #elif defined MULTICOMP && defined FD
                        for(size_t a = 0; a < nSpecies; ++a)
                        {
                            species[a][i][j][k].rho = mixture[i][j][k].rho * gas->massFraction(gas->speciesIndex(speciesName[a]));
                        }
                    #endif
                }
            }
        }
    }

    fill_BC();
}

void LBM::Collide(){    
     #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){    
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E){   
                    #ifndef MULTICOMP
                        double cv = gas_const / (gamma - 1.0);
                        double cp = cv + gas_const;
                        double theta = gas_const*mixture[i][j][k].temp;
                        double mu = nu*mixture[i][j][k].rho;
                        double conduc_coeff = mu*cp/prtl;
                        // std::cout << "alpha : " << conduc_coeff/(cp * mixture[i][j][k].rho) << std::endl;
                        double enthalpy = cp * mixture[i][j][k].temp; // H = Cp * T = (Cv + 1) * T
                    
                    #else
                        // Create Cantera object
                        int rank = omp_get_thread_num();
                        auto gas = sols[rank]->thermo();
                        std::vector <double> Y (gas->nSpecies());
                        for(size_t a = 0; a < nSpecies; ++a) Y [gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                        gas->setMassFractions(&Y[0]);
                        gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));

                        // Get macropscopic properties
                        double cp = units.cp(gas->cp_mass());
                        double enthalpy = units.energy_mass(gas->enthalpy_mass());
                        std::vector <double> part_enthalpy(gas->nSpecies());
                        gas->getPartialMolarEnthalpies(&part_enthalpy[0]);

                        auto trans = sols[rank]->transport();
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

                    for(int p=0; p < 3; ++p){
                        for(int q=0; q < 3; ++q){
                            eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 
                            eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + mixture[i][j][k].p*velocity[p]*velocity[q];
                        }  
                    }

                    // double delYx[nSpecies], delYy[nSpecies], delYz[nSpecies];
                    // for(size_t a = 0; a < nSpecies; ++a){
                    //     delYx[a] = ((species[a][i+1][j][k].rho / mixture[i+1][j][k].rho) - (species[a][i-1][j][k].rho / mixture[i-1][j][k].rho)) / (dx);
                    //     delYy[a] = ((species[a][i][j+1][k].rho / mixture[i][j+1][k].rho) - (species[a][i][j-1][k].rho / mixture[i][j-1][k].rho)) / (dy);
                    //     delYz[a] = ((species[a][i][j][k+1].rho / mixture[i][j][k+1].rho) - (species[a][i][j][k-1].rho / mixture[i][j][k-1].rho)) / (dz);
                    // }

                    double delQdevx, delQdevy, delQdevz;
                    delQdevx = mixture[i][j][k].dQdevx ;//FD_limiterVanleer( mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1-3*gas_const*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1-3*gas_const*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u), dx) ;
                    delQdevy = mixture[i][j][k].dQdevy ;//FD_limiterVanleer( mixture[i][j+1][k].rho*mixture[i][j+1][k].u*(1-3*gas_const*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i][j-1][k].rho*mixture[i][j-1][k].u*(1-3*gas_const*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].u), dy) ;
                    delQdevz = mixture[i][j][k].dQdevz ;//FD_limiterVanleer( mixture[i][j][k+1].rho*mixture[i][j][k+1].u*(1-3*gas_const*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i][j][k-1].rho*mixture[i][j][k-1].u*(1-3*gas_const*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].u), dz) ;

                    double q_diff [3] = {0.0, 0.0, 0.0};
                    double q_corr [3] = {0.0, 0.0, 0.0};
                    #ifdef MULTICOMP
                    for (size_t a = 0; a < nSpecies; ++a){
                        int speciesIdx = gas->speciesIndex(speciesName[a]);
                        double mmass = gas->molecularWeight(speciesIdx);
                        q_diff[0] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].Vdiff_x);
                        q_diff[1] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].Vdiff_y);
                        q_diff[2] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].Vdiff_z);
                    
                        q_corr[0] += 0.5 * (omega1-2.0)/(omega1-omega) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * species[a][i][j][k].delYx;
                        q_corr[1] += 0.5 * (omega1-2.0)/(omega1-omega) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * species[a][i][j][k].delYy;
                        q_corr[2] += 0.5 * (omega1-2.0)/(omega1-omega) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * species[a][i][j][k].delYz;              
                    }
                    #endif

                    double str_heat_flux[3] ={  mixture[i][j][k].energy_flux[0] - velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) - velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) - velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) - 0.5*dt_sim*velocity[0]*delQdevx + q_diff[0] + q_corr[0],   //
                                                mixture[i][j][k].energy_flux[1] - velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) - velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) - velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) - 0.5*dt_sim*velocity[1]*delQdevy + q_diff[1] + q_corr[1],   //
                                                mixture[i][j][k].energy_flux[2] - velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) - velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) - velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) - 0.5*dt_sim*velocity[2]*delQdevz + q_diff[2] + q_corr[2]} ; //

                    double corr[3] = {  dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevx,
                                        dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevy,
                                        dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevz};

                    // double corr[3] = {0.0, 0.0, 0.0};

                    for (int l = 0; l < npop; ++l){
                        // ------------- Mass and Momentum collision -----------------------------
                        double feq = calculate_feq(l, mixture[i][j][k].rho, velocity, theta, corr);
                        double geq = calculate_geq(l, mixture[i][j][k].rhoe, eq_heat_flux, eq_R_tensor, theta);
                        double gstr = calculate_geq(l, mixture[i][j][k].rhoe, str_heat_flux, eq_R_tensor, theta);

                        mixture[i][j][k].fpc[l] = (1.0-omega)*mixture[i][j][k].f[l] + omega*feq;
                        mixture[i][j][k].gpc[l] = mixture[i][j][k].g[l] + omega1*(geq-mixture[i][j][k].g[l]) + (omega-omega1)*(gstr-mixture[i][j][k].g[l]);
                    }     
                }
            }
        }
    }
}


void LBM::Streaming(){
    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            for(int k = 0; k<Nz; ++k){
                if(mixture[i][j][k].type==TYPE_F){
                    for (int l=0; l < npop; ++l){
                        int i_nb, j_nb, k_nb;
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        //---- Solid Boundary Condition ----------------------
                        if(mixture[i_nb][j_nb][k_nb].type==TYPE_S){
                            mixture[i][j][k].f[l] = mixture[i][j][k].fpc[opposite[l]];
                            // mixture[i][j][k].g[l] = mixture[i][j][k].gpc[opposite[l]];
                            #if defined MULTICOMP && !defined FD
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = species[a][i][j][k].fpc[opposite[l]];
                            #endif
                        }
                        //---- Adiabatic Wall --------------------------
                        else if (mixture[i_nb][j_nb][k_nb].type==TYPE_A){
                            mixture[i][j][k].f[l] = mixture[i][j][k].fpc[opposite[l]];
                            mixture[i][j][k].g[l] = mixture[i][j][k].gpc[opposite[l]];
                            #if defined MULTICOMP && !defined FD
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = species[a][i][j][k].fpc[opposite[l]];
                            #endif

                        }
                        //---- Inlet/Outlet Boundary Condition ---------------
                        else if (mixture[i_nb][j_nb][k_nb].type==TYPE_E){
                            mixture[i][j][k].f[l] = mixture[i_nb][j_nb][k_nb].fpc[l];
                            mixture[i][j][k].g[l] = mixture[i_nb][j_nb][k_nb].gpc[l];
                        }
                        //---- Periodic Boundary Condition --------------------
                        else {
                            /* Alternative Periodic Code
                            if (i_nb < 1) i_nb = Nx-2;
                            else if(i_nb > Nx-2) i_nb = 1;

                            if (j_nb < 1) j_nb = Ny-2;
                            else if(j_nb > Ny-2) j_nb = 1;

                            if (k_nb < 1) k_nb = Nz-2;
                            else if(k_nb > Nz-2) k_nb = 1;
                            */

                            
                            i_nb = ((i_nb - 1 + (Nx-2)) % (Nx-2)) + 1;
                            j_nb = ((j_nb - 1 + (Ny-2)) % (Ny-2)) + 1;
                            k_nb = ((k_nb - 1 + (Nz-2)) % (Nz-2)) + 1;
                            mixture[i][j][k].f[l] = mixture[i_nb][j_nb][k_nb].fpc[l];
                            mixture[i][j][k].g[l] = mixture[i_nb][j_nb][k_nb].gpc[l];
                            #if defined MULTICOMP && !defined FD
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = species[a][i_nb][j_nb][k_nb].fpc[l];
                            #endif
                            
                        }
                    }
                }
            }
        }
    }

    // TMS Boundary Conditions
    #ifdef PARALLEL 
    #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            for(int k = 0; k<Nz; ++k){
                if(mixture[i][j][k].type==TYPE_F){
                    int i_nb, j_nb, k_nb;
                    double f_tgt[npop];
                    double g_tgt[npop];
                    double f_loc[npop];
                    double g_loc[npop];
                    bool interface_node = false; 
                    #ifdef MULTICOMP
                        double fa_tgt[nSpecies][npop];
                        double fa_loc[nSpecies][npop];
                    #endif

                    // check interface node                         
                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if(mixture[i_nb][j_nb][k_nb].type==TYPE_S){    
                            interface_node = true;
                            break;
                        }                            
                    }

                    // step 1: calculate f_tgt, g_tft, and fa_tgt
                    if(interface_node == true){
                        double rho_bb = 0.0;
                        #ifdef MULTICOMP
                        double Y_bb[nSpecies];
                        #endif

                        for (int l=0; l < npop; ++l){
                            rho_bb += mixture[i][j][k].f[l];
                            
                        #ifdef MULTICOMP
                            for(size_t a = 0; a < nSpecies; ++a)
                                Y_bb[a] += species[a][i][j][k].f[l] / rho_bb;
                        #endif
                        }

                        double vel_tgt[3] = {mixture[(int) i_nb][(int) j_nb][(int) k_nb].u,
                                            mixture[(int) i_nb][(int) j_nb][(int) k_nb].v,
                                            mixture[(int) i_nb][(int) j_nb][(int) k_nb].w};
                        double T_tgt = mixture[(int) i_nb][(int) j_nb][(int) k_nb].temp;                      
                        
                        #ifndef MULTICOMP
                            double p_tgt = rho_bb * gas_const * T_tgt;
                            double theta = gas_const*T_tgt;
                            double cv = gas_const / (gamma - 1.0);
                            double cp = cv + gas_const;
                            double internal_energy = cv * T_tgt;
                            double enthalpy = cp * T_tgt; // H = Cp * T = (Cv + 1) * T
                        #else
                            int rank = omp_get_thread_num();
                            auto gas = sols[rank]->thermo();
                            std::vector <double> Y (gas->nSpecies());
                            for(size_t a = 0; a < nSpecies; ++a)
                                Y[gas->speciesIndex(speciesName[a])] = Y_bb[a];
                            
                            gas->setMassFractions(&Y[0]);
                            gas->setState_TD(units.si_temp(T_tgt), units.si_rho(rho_bb));

                            double internal_energy = units.energy_mass(gas->intEnergy_mass());
                            double enthalpy = units.energy_mass(gas->enthalpy_mass());

                            double p_tgt = units.p(gas->pressure());
                            double theta = units.energy_mass(gas->RT()/gas->meanMolecularWeight());      
                        #endif

                        double rhoe = rho_bb*(internal_energy + 0.5 * v_sqr(vel_tgt[0], vel_tgt[1], vel_tgt[2]));
                        double total_enthalpy = enthalpy + 0.5 * v_sqr(vel_tgt[0], vel_tgt[1], vel_tgt[2]);
                        double eq_heat_flux[3] = {  total_enthalpy*rho_bb*vel_tgt[0] ,
                                                    total_enthalpy*rho_bb*vel_tgt[1] ,
                                                    total_enthalpy*rho_bb*vel_tgt[2] };
                        double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                    {0., 0., 0.},
                                                    {0., 0., 0.}};
                        double eq_R_tensor[3][3] = {{0., 0., 0.},    // second-order moment of g
                                                    {0., 0., 0.},
                                                    {0., 0., 0.}};

                        for(int p=0; p < 3; ++p){
                            for(int q=0; q < 3; ++q){
                                eq_p_tensor[p][q] = (p==q) ? p_tgt+rho_bb*vel_tgt[p]*vel_tgt[q] : rho_bb*vel_tgt[p]*vel_tgt[q]; 
                                eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + p_tgt*vel_tgt[p]*vel_tgt[q];
                            }  
                        }
                        
                        double corr[3]= {0, 0, 0};

                        for (int l=0; l < npop; ++l){   
                            f_tgt[l] = calculate_feq(l, rho_bb, vel_tgt, theta, corr);
                            g_tgt[l] = calculate_geq(l, rhoe, eq_heat_flux, eq_R_tensor, theta);

                            #ifdef MULTICOMP
                            for(size_t a = 0; a < nSpecies; ++a){
                                theta = units.energy_mass(gas->RT() / gas->molecularWeight(gas->speciesIndex(speciesName[a])) );
                                fa_tgt[a][l] = calculate_feq(l, rho_bb*Y_bb[a], vel_tgt, theta, corr);
                            }
                            #endif
                        }
                    }
                    
                    // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
                    if(interface_node == true) {
                        double rho_loc = 0.0;
                        double rhou_loc = 0.0;
                        double rhov_loc = 0.0;
                        double rhow_loc = 0.0;
                        double rhoe_loc = 0.0;
                        #ifdef MULTICOMP
                            double rhoa_loc[nSpecies] ;
                            double rhoua_loc[nSpecies] ;
                            double rhova_loc[nSpecies] ;
                            double rhowa_loc[nSpecies] ;
                        #endif

                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
                                rho_loc += f_tgt[l];
                                rhou_loc += f_tgt[l]*cx[l];
                                rhov_loc += f_tgt[l]*cy[l];
                                rhow_loc += f_tgt[l]*cz[l];
                                rhoe_loc += g_tgt[l];

                                #ifdef MULTICOMP
                                for(size_t a=0; a < nSpecies; ++a){
                                    rhoa_loc[a] += fa_tgt[a][l];
                                    rhoua_loc[a] += fa_tgt[a][l]*cx[l];
                                    rhova_loc[a] += fa_tgt[a][l]*cy[l];
                                    rhowa_loc[a] += fa_tgt[a][l]*cz[l];
                                }
                                #endif
                            }
                            else{
                                rho_loc += mixture[i][j][k].f[l];
                                rhou_loc += mixture[i][j][k].f[l]*cx[l];
                                rhov_loc += mixture[i][j][k].f[l]*cy[l];
                                rhow_loc += mixture[i][j][k].f[l]*cz[l];
                                rhoe_loc += mixture[i][j][k].g[l];

                                #ifdef MULTICOMP
                                    for(size_t a=0; a < nSpecies; ++a){
                                        rhoa_loc[a] += species[a][i][j][k].f[l];
                                        rhoua_loc[a] += species[a][i][j][k].f[l]*cx[l];
                                        rhova_loc[a] += species[a][i][j][k].f[l]*cy[l];
                                        rhowa_loc[a] += species[a][i][j][k].f[l]*cz[l];
                                    }
                                #endif
                            }
                        }

                        double vel_loc[3];
                        vel_loc[0] = rhou_loc / rho_loc;
                        vel_loc[1] = rhov_loc / rho_loc;
                        vel_loc[2] = rhow_loc / rho_loc;

                        #ifdef MULTICOMP
                            double vela_loc[nSpecies][3];
                            for(size_t a=0; a < nSpecies; ++a){
                                vela_loc[a][0] = rhoua_loc[a] / rhoa_loc[a];
                                vela_loc[a][1] = rhova_loc[a] / rhoa_loc[a];
                                vela_loc[a][2] = rhowa_loc[a] / rhoa_loc[a];
                            }
                        #endif
                        
                        double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                        
                        #ifndef MULTICOMP
                            double cv = gas_const / (gamma - 1.0);
                            double T_loc = internalEnergy / cv;                        
                            
                            double p_loc = rho_loc * gas_const * T_loc;
                            double theta = gas_const*T_loc;
                            double cp = cv + gas_const;
                            double internal_energy = cv * T_loc;
                            double enthalpy = cp * T_loc; // H = Cp * T = (Cv + 1) * T
                        #else
                            int rank = omp_get_thread_num();
                            auto gas = sols[rank]->thermo();
                            std::vector <double> Y_loc (gas->nSpecies());

                            for(size_t a = 0; a < nSpecies; ++a){
                                Y_loc[gas->speciesIndex(speciesName[a])] = rhoa_loc[a] / rho_loc;
                                std::cout << rhoa_loc[a] / rho_loc << std::endl;
                            }
                            gas->setMassFractions(&Y_loc[0]);
                            gas->setState_UV(units.si_energy_mass(internalEnergy), 1.0 / units.si_rho(rho_loc));

                            double internal_energy = units.energy_mass(gas->intEnergy_mass());
                            double enthalpy = units.energy_mass(gas->enthalpy_mass());

                            double p_loc = units.p(gas->pressure());
                            double theta = units.energy_mass(gas->RT()/gas->meanMolecularWeight());
                        #endif

                        double rhoe = rho_loc*(internal_energy + 0.5 * v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]));
                        double total_enthalpy = enthalpy + 0.5 * v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                        double eq_heat_flux[3] = {  total_enthalpy*rho_loc*vel_loc[0],
                                                    total_enthalpy*rho_loc*vel_loc[1],
                                                    total_enthalpy*rho_loc*vel_loc[2]};
                        double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                    {0., 0., 0.},
                                                    {0., 0., 0.}};
                        double eq_R_tensor[3][3] = {{0., 0., 0.},    // second-order moment of g
                                                    {0., 0., 0.},
                                                    {0., 0., 0.}};

                        for(int p=0; p < 3; ++p){
                            for(int q=0; q < 3; ++q){
                                eq_p_tensor[p][q] = (p==q) ? p_loc+rho_loc*vel_loc[p]*vel_loc[q] : rho_loc*vel_loc[p]*vel_loc[q]; 
                                eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + p_loc*vel_loc[p]*vel_loc[q];
                            }  
                        }
                        
                        double corr[3]= {0, 0, 0};
                        
                        for (int l=0; l < npop; ++l){
                            f_loc[l] = calculate_feq(l, rho_loc, vel_loc, theta, corr);
                            g_loc[l] = calculate_geq(l, rhoe, eq_heat_flux, eq_R_tensor, theta);

                            // #ifdef MULTICOMP
                            // for(int a = 0; a < nSpecies; ++a)
                            // {
                            //     theta = units.energy_mass(gas->RT() / gas->molecularWeight(gas->speciesIndex(speciesName[a])) );
                            //     fa_loc[a][l] = calculate_feq(l, rhoa_loc[a], vela_loc[a], theta, corr);
                            // }
                            // #endif
                        }

                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
                                mixture[i][j][k].f[l] = 2*f_tgt[l] - f_loc[l];
                                mixture[i][j][k].g[l] = 2*g_tgt[l] - g_loc[l];

                                // #ifdef MULTICOMP
                                // for(int a = 0; a < nSpecies; ++a)
                                //     species[a][i][j][k].f[l] = 2*fa_tgt[a][l] - fa_loc[a][l];
                                // #endif
                            }
                            else{
                                mixture[i][j][k].f[l] = f_tgt[l] + mixture[i][j][k].f[l] - f_loc[l];
                                mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];

                                // #ifdef MULTICOMP
                                // for(int a = 0; a < nSpecies; ++a)
                                //     species[a][i][j][k].f[l] = fa_tgt[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];
                                // #endif
                            }
                        }
                    }

                } 
            }
        }
    }
}

#ifdef MULTICOMP
void LBM::Collide_Species(){
    #ifdef PARALLEL 
    #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){    
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E){   
                    double mass_frac[nSpecies];         
                    double tau_ab[nSpecies][nSpecies];
                    double D_ab[nSpecies][nSpecies];
                    double invtau_a[nSpecies];
                    double lambda[nSpecies];
                    double beta[nSpecies];
                    double mmass[nSpecies]; 
                    Eigen::SparseMatrix<double> matrix_A(nSpecies, nSpecies);
                    Eigen::VectorXd vector_u(nSpecies);
                    Eigen::VectorXd vector_v(nSpecies);
                    Eigen::VectorXd vector_w(nSpecies);
                    Eigen::VectorXd du(nSpecies);
                    Eigen::VectorXd dv(nSpecies);
                    Eigen::VectorXd dw(nSpecies);

                    int rank = omp_get_thread_num();
                    auto gas = sols[rank]->thermo();               
                    std::vector<double> Y (gas->nSpecies());
                    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));
                    
                    for(size_t a = 0; a < nSpecies; ++a) {
                        mass_frac[a] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                        mmass[a] = gas->molecularWeight(gas->speciesIndex(speciesName[a]));
                    }

                    double mix_mmass = gas->meanMolecularWeight();
                    double theta = gas->RT(); // universal_gas_const * temp

                    int ld = gas->nSpecies();
                    std::vector<double> d(ld * ld);
                    auto trans = sols[rank]->transport();
                    trans->getBinaryDiffCoeffs(ld, &d[0]);

                    for(size_t a = 0; a < nSpecies; ++a){                       
                        for(size_t b = 0; b < nSpecies; ++b){
                            // get diffusion coefficient
                            D_ab[a][b] =   d[ld*gas->speciesIndex(speciesName[b]) + gas->speciesIndex(speciesName[a])] ;
                            // std::cout << "D_ab " << gas->speciesName(gas->speciesIndex(speciesName[a])) << "-" << gas->speciesName(gas->speciesIndex(speciesName[b])) << " : " << D_ab[a][b] << std::endl;
                            
                            tau_ab[a][b] = units.t( (mmass[a]*mmass[b]/(mix_mmass*theta)) * D_ab[a][b] );
                            // std::cout << " tau_ " <<  gas->speciesName(gas->speciesIndex(speciesName[a])) << "-" << gas->speciesName(gas->speciesIndex(speciesName[b])) << " : " << tau_ab[a][b] << std::endl;
                        }
                    }
                    
                    for(size_t a = 0; a < nSpecies; ++a){    
                        invtau_a[a] = 0.0;
                        for(size_t b = 0; b < nSpecies; ++b){                                
                            invtau_a[a] += mass_frac[b] / tau_ab[a][b];  
                        }
                    }

                    if(diffModel == 0){  // Stefan-Maxwell Diffusion Model
                        for(size_t a = 0; a < nSpecies; ++a){
                            for(size_t b = 0; b < nSpecies; ++b){
                                if (a == b) matrix_A.insert(a, b) = 1.0 + dt_sim * invtau_a[a] / 2.0 - dt_sim/2.0*mass_frac[a]/tau_ab[a][b];
                                else matrix_A.insert(a, b) = - dt_sim / 2.0 * mass_frac[b] / tau_ab[a][b];
                            }

                            vector_u(a) = species[a][i][j][k].u - mixture[i][j][k].u;
                            vector_v(a) = species[a][i][j][k].v - mixture[i][j][k].v;
                            vector_w(a) = species[a][i][j][k].w - mixture[i][j][k].w;

                            beta[a] = dt_sim / (2*(1/invtau_a[a]) + dt_sim);
                        }  
                        
                        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
                        solver.analyzePattern(matrix_A);
                        solver.factorize(matrix_A);                
                        solver.compute(matrix_A);

                        du = solver.solve(vector_u);
                        dv = solver.solve(vector_v);
                        dw = solver.solve(vector_w);

                        for(size_t a = 0; a < nSpecies; ++a){
                            species[a][i][j][k].Vdiff_x = du[a];
                            species[a][i][j][k].Vdiff_y = dv[a];
                            species[a][i][j][k].Vdiff_z = dw[a];
                        }
                    }
                    else{   // Mixture-averaged Diffusion Model
                        double denom = 0.0;
                        for (size_t b = 0; b < nSpecies; ++b)
                            denom += mass_frac[b] * invtau_a[b];

                        for (size_t a = 0; a < nSpecies; ++a){                                
                            lambda[a] = mass_frac[a] * invtau_a[a] / denom;
                            beta[a] = dt_sim / (2*(1/invtau_a[a]) + dt_sim);
                        }

                        double betalamda = {0.0};
                        double betaYu[nSpecies][ndim] = {0.0};
                        for (size_t a = 0; a < nSpecies; ++a){
                            if (mass_frac[a] == 0.0) continue;

                            betalamda += beta[a] * lambda[a];
                            for (size_t b = 0; b < nSpecies; ++b){
                                betaYu[a][0] += beta[b] * mass_frac[b] / mass_frac[a] * (species[b][i][j][k].u - mixture[i][j][k].u);
                                betaYu[a][1] += beta[b] * mass_frac[b] / mass_frac[a] * (species[b][i][j][k].v - mixture[i][j][k].v);
                                betaYu[a][2] += beta[b] * mass_frac[b] / mass_frac[a] * (species[b][i][j][k].w - mixture[i][j][k].w);
                            }
                        }

                        for (size_t a = 0; a < nSpecies; ++a){
                            species[a][i][j][k].Vdiff_x = (1.0 - beta[a]) * (species[a][i][j][k].u - mixture[i][j][k].u) + (1.0 - beta[a]) * lambda[a] / (1.0 - betalamda) * betaYu[a][0] ;
                            species[a][i][j][k].Vdiff_y = (1.0 - beta[a]) * (species[a][i][j][k].v - mixture[i][j][k].v) + (1.0 - beta[a]) * lambda[a] / (1.0 - betalamda) * betaYu[a][1] ;
                            species[a][i][j][k].Vdiff_z = (1.0 - beta[a]) * (species[a][i][j][k].w - mixture[i][j][k].w) + (1.0 - beta[a]) * lambda[a] / (1.0 - betalamda) * betaYu[a][2] ;
                        }
                    }

                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};

                    // std::vector <double> w_dot(gas->nSpecies());    // mole density rate [kmol/m3/s]
                    // auto kinetics = sols[rank]->kinetics();
// 
                    // size_t p = 0;
                    // double R_a[gas->nSpecies()] = {0.0};
                    // double rho_a[gas->nSpecies()] = {0.0};
                    // do{
                        // kinetics->getNetProductionRates(&w_dot[0]); 
                        // for (int a = 0; a < (int) gas->nSpecies(); ++a){                            
                            // if (w_dot[a] == 0.0)  // Check rate == 0, no reaction take place
                                // continue;
// 
                            // #pragma omp critical
                            // { 
                                // bool new_species = true;
                                // for (size_t b = 0; b < nSpecies; ++b){
                                    // if(gas->speciesName(a) == speciesName[b]){
                                        // new_species = false;
                                        // break;
                                    // }
                                // }
// 
                                // if (new_species){   // adding product to the LBM variables
                                    // speciesName.push_back(gas->speciesName(a));
                                    // nSpecies++;
// 
                                    // allocate memory for new species
                                    // this->species.resize(nSpecies);
                                    // this->species[nSpecies-1] = new SPECIES **[this->Nx];
                                    // for (int i = 0; i < this->Nx; ++i){
                                        // this->species[nSpecies-1][i] = new SPECIES *[this->Ny];
                                        // for (int j = 0; j < this->Ny; ++j){
                                            // this->species[nSpecies-1][i][j] = new SPECIES [this->Nz];
                                        // }
                                    // }                            
                                // } 
                            // }
// 
                            // size_t idx_species = 0;
                            // for (size_t b = 0; b < nSpecies; ++b){
                                // if(gas->speciesName(a) == speciesName[b]){
                                    // idx_species = b;
                                    // break;
                                // }
                            // }
                            // 
                            // double rho_dot = units.rho_dot(w_dot[a] * gas->molecularWeight(a));
                            // R_a[idx_species] += rho_dot;
                            // rho_a[idx_species] += dt_sim/2.0 * rho_dot;
                        // }
// 
                        // for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = (species[a][i][j][k].rho+rho_a[a]) / mixture[i][j][k].rho;
                        // gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));
                        // gas->setMassFractions(&Y[0]);
                        // p++;
                    // }while(p < 2);
                    
                    // for(size_t a = 0; a < nSpecies; ++a) 
                    //     R_a[a] = dt_sim/2.0 * R_a[a];

                    theta = units.energy_mass(gas->RT()/gas->meanMolecularWeight());
                    
                    double corr[nSpecies][3] = {0};

                    // for(size_t a = 0; a < nSpecies; ++a){
                    //     double delQdevx = FD_limiterVanleer( species[a][i+1][j][k].rho*species[a][i+1][j][k].u*(1-3*gas_const*mixture[i+1][j][k].temp)-species[a][i+1][j][k].rho*cb(species[a][i+1][j][k].u), species[a][i][j][k].rho*species[a][i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-species[a][i][j][k].rho*cb(species[a][i][j][k].u), species[a][i-1][j][k].rho*species[a][i-1][j][k].u*(1-3*gas_const*mixture[i-1][j][k].temp)-species[a][i-1][j][k].rho*cb(species[a][i-1][j][k].u), dx) ;
                    //     double delQdevy = FD_limiterVanleer( species[a][i][j+1][k].rho*species[a][i][j+1][k].u*(1-3*gas_const*mixture[i][j+1][k].temp)-species[a][i][j+1][k].rho*cb(species[a][i][j+1][k].u), species[a][i][j][k].rho*species[a][i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-species[a][i][j][k].rho*cb(species[a][i][j][k].u), species[a][i][j-1][k].rho*species[a][i][j-1][k].u*(1-3*gas_const*mixture[i][j-1][k].temp)-species[a][i][j-1][k].rho*cb(species[a][i][j-1][k].u), dy) ;
                    //     double delQdevz = FD_limiterVanleer( species[a][i][j][k+1].rho*species[a][i][j][k+1].u*(1-3*gas_const*mixture[i][j][k+1].temp)-species[a][i][j][k+1].rho*cb(species[a][i][j][k+1].u), species[a][i][j][k].rho*species[a][i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-species[a][i][j][k].rho*cb(species[a][i][j][k].u), species[a][i][j][k-1].rho*species[a][i][j][k-1].u*(1-3*gas_const*mixture[i][j][k-1].temp)-species[a][i][j][k-1].rho*cb(species[a][i][j][k-1].u), dz) ;
                    
                    //     corr[a][0] = dt_sim*(2.0-2.0*beta[a])/(2.0*species[a][i][j][k].rho*2.0*beta[a])*delQdevx;
                    //     corr[a][1] = dt_sim*(2.0-2.0*beta[a])/(2.0*species[a][i][j][k].rho*2.0*beta[a])*delQdevy;
                    //     corr[a][2] = dt_sim*(2.0-2.0*beta[a])/(2.0*species[a][i][j][k].rho*2.0*beta[a])*delQdevz;
                    // }

                    

                    for (int l = 0; l < npop; ++l){
                        double feq[nSpecies], fstr[nSpecies];

                        // double feqM = calculate_feq(l, 1.0, velocity, theta, corr);

                        for (size_t a = 0; a < nSpecies; ++a){
                            double velocity_spec[3] = { velocity[0] + species[a][i][j][k].Vdiff_x,
                                                        velocity[1] + species[a][i][j][k].Vdiff_y, 
                                                        velocity[2] + species[a][i][j][k].Vdiff_z};

                            theta = units.energy_mass( gas->RT()/gas->molecularWeight(gas->speciesIndex(speciesName[a])) );
                            feq[a] = calculate_feq(l, species[a][i][j][k].rho, velocity, theta, corr[a]);
                            fstr[a] = calculate_feq(l, species[a][i][j][k].rho, velocity_spec, theta, corr[a]);
                            // freact[a] = calculate_feq(l, species[a][i][j][k].rho_dot, velocity, theta, corr);
                        }

                        for (size_t a = 0; a < nSpecies; ++a){
                            double F_a = 0.0;

                            if (diffModel == 0) // Stefan-Maxwell Diffusion
                                for(size_t b = 0; b < nSpecies; ++b)
                                    F_a += mass_frac[a]/tau_ab[a][b]*(feq[b]-fstr[b]);
                            else    // Mixture-averaged Diffusion
                                for(size_t b = 0; b < nSpecies; ++b)
                                    F_a += lambda[a]*invtau_a[b]*(feq[b]-fstr[b]);                                

                            // species[a][i][j][k].fpc[l] = species[a][i][j][k].f[l] + 2.0*beta[a]*(feq[a]-species[a][i][j][k].f[l]) + dt_sim*(beta[a]-1.0)*F_a + dt_sim*freact[a];
                            // species[a][i][j][k].fpc[l] = species[a][i][j][k].f[l] + 2.0*beta[a]*(feq[a]-species[a][i][j][k].f[l]) + dt_sim*(beta[a]-1.0)*F_a + feqM*R_a[a];
                            species[a][i][j][k].fpc[l] = species[a][i][j][k].f[l] + 2.0*beta[a]*(feq[a]-species[a][i][j][k].f[l]) + dt_sim*(beta[a]-1.0)*F_a;
                        }
                    }   
                }
            }
        }
    }
}

#endif

void LBM::fill_BC()
{
    for(int i = 0; i < Nx; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){
                if(mixture[i][j][k].type == TYPE_F){
                    int i_nb = i - 1;
                    int i_pb = i + 1;
                    int j_nb = j - 1;
                    int j_pb = j + 1;
                    int k_nb = k - 1;
                    int k_pb = k + 1;

                    // Adiabatic Boundary Condition
                    if(mixture[i_nb][j][k].type == TYPE_A){
                        mixture[i_nb][j][k].rho = mixture[i][j][k].rho;
                        mixture[i_nb][j][k].u   = mixture[i][j][k].u;
                        mixture[i_nb][j][k].v   = mixture[i][j][k].v;
                        mixture[i_nb][j][k].w   = mixture[i][j][k].w;
                        mixture[i_nb][j][k].temp= mixture[i][j][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_nb][j][k].rho = species[a][i][j][k].rho;
                    }
                    if(mixture[i_pb][j][k].type == TYPE_A){
                        mixture[i_pb][j][k].rho = mixture[i][j][k].rho;
                        mixture[i_pb][j][k].u   = mixture[i][j][k].u;
                        mixture[i_pb][j][k].v   = mixture[i][j][k].v;
                        mixture[i_pb][j][k].w   = mixture[i][j][k].w;
                        mixture[i_pb][j][k].temp= mixture[i][j][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_pb][j][k].rho = species[a][i][j][k].rho;
                    }
                    if(mixture[i][j_nb][k].type == TYPE_A){
                        mixture[i][j_nb][k].rho = mixture[i][j][k].rho;
                        mixture[i][j_nb][k].u   = mixture[i][j][k].u;
                        mixture[i][j_nb][k].v   = mixture[i][j][k].v;
                        mixture[i][j_nb][k].w   = mixture[i][j][k].w;
                        mixture[i][j_nb][k].temp= mixture[i][j][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_nb][k].rho = species[a][i][j][k].rho;
                    }
                    if(mixture[i][j_pb][k].type == TYPE_A){
                        mixture[i][j_pb][k].rho = mixture[i][j][k].rho;
                        mixture[i][j_pb][k].u   = mixture[i][j][k].u;
                        mixture[i][j_pb][k].v   = mixture[i][j][k].v;
                        mixture[i][j_pb][k].w   = mixture[i][j][k].w;
                        mixture[i][j_pb][k].temp= mixture[i][j][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_pb][k].rho = species[a][i][j][k].rho;
                    }
                    if(mixture[i][j][k_nb].type == TYPE_A){
                        mixture[i][j][k_nb].rho = mixture[i][j][k].rho;
                        mixture[i][j][k_nb].u   = mixture[i][j][k].u;
                        mixture[i][j][k_nb].v   = mixture[i][j][k].v;
                        mixture[i][j][k_nb].w   = mixture[i][j][k].w;
                        mixture[i][j][k_nb].temp= mixture[i][j][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_nb].rho = species[a][i][j][k].rho;
                    }
                    if(mixture[i][j][k_pb].type == TYPE_A){
                        mixture[i][j][k_pb].rho = mixture[i][j][k].rho;
                        mixture[i][j][k_pb].u   = mixture[i][j][k].u;
                        mixture[i][j][k_pb].v   = mixture[i][j][k].v;
                        mixture[i][j][k_pb].w   = mixture[i][j][k].w;
                        mixture[i][j][k_pb].temp= mixture[i][j][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_pb].rho = species[a][i][j][k].rho;
                    }
                    
                    // Periodic Boundary Condition 
                    if(mixture[i_nb][j][k].type == TYPE_P){
                        mixture[i_nb][j][k].rho = mixture[Nx-2][j][k].rho;
                        mixture[i_nb][j][k].u   = mixture[Nx-2][j][k].u;
                        mixture[i_nb][j][k].v   = mixture[Nx-2][j][k].v;
                        mixture[i_nb][j][k].w   = mixture[Nx-2][j][k].w;
                        mixture[i_nb][j][k].temp= mixture[Nx-2][j][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_nb][j][k].rho = species[a][Nx-2][j][k].rho;
                    }
                    if(mixture[i_pb][j][k].type == TYPE_P){
                        mixture[i_pb][j][k].rho = mixture[1][j][k].rho;
                        mixture[i_pb][j][k].u   = mixture[1][j][k].u;
                        mixture[i_pb][j][k].v   = mixture[1][j][k].v;
                        mixture[i_pb][j][k].w   = mixture[1][j][k].w;
                        mixture[i_pb][j][k].temp= mixture[1][j][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_pb][j][k].rho = species[a][1][j][k].rho;
                    }
                    if(mixture[i][j_nb][k].type == TYPE_P){
                        mixture[i][j_nb][k].rho = mixture[i][Ny-2][k].rho;
                        mixture[i][j_nb][k].u   = mixture[i][Ny-2][k].u;
                        mixture[i][j_nb][k].v   = mixture[i][Ny-2][k].v;
                        mixture[i][j_nb][k].w   = mixture[i][Ny-2][k].w;
                        mixture[i][j_nb][k].temp= mixture[i][Ny-2][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_nb][k].rho = species[a][i][Ny-2][k].rho;
                    }
                    if(mixture[i][j_pb][k].type == TYPE_P){
                        mixture[i][j_pb][k].rho = mixture[i][1][k].rho;
                        mixture[i][j_pb][k].u   = mixture[i][1][k].u;
                        mixture[i][j_pb][k].v   = mixture[i][1][k].v;
                        mixture[i][j_pb][k].w   = mixture[i][1][k].w;
                        mixture[i][j_pb][k].temp= mixture[i][1][k].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_pb][k].rho = species[a][i][1][k].rho;
                    }
                    if(mixture[i][j][k_nb].type == TYPE_P){
                        mixture[i][j][k_nb].rho = mixture[i][j][Nz-2].rho;
                        mixture[i][j][k_nb].u   = mixture[i][j][Nz-2].u;
                        mixture[i][j][k_nb].v   = mixture[i][j][Nz-2].v;
                        mixture[i][j][k_nb].w   = mixture[i][j][Nz-2].w;
                        mixture[i][j][k_nb].temp= mixture[i][j][Nz-2].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_nb].rho = species[a][i][j][Nz-2].rho;
                    }
                    if(mixture[i][j][k_pb].type == TYPE_P){
                        mixture[i][j][k_pb].rho = mixture[i][j][1].rho;
                        mixture[i][j][k_pb].u   = mixture[i][j][1].u;
                        mixture[i][j][k_pb].v   = mixture[i][j][1].v;
                        mixture[i][j][k_pb].w   = mixture[i][j][1].w;
                        mixture[i][j][k_pb].temp= mixture[i][j][1].temp;
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_pb].rho = species[a][i][j][1].rho;
                    }


                }
            }
        }
    }
}

void LBM::run(int nstep, int tout)
{ 
    this->nstep = nstep;
    this->tout = tout;

    std::cout << "  Setup Done" << std::endl;

    // initialize the distribution function 
    std::cout << "  Initialization ..." << std::endl;
    Init();  
    std::cout << "  Initialization Done" << std::endl;
    
    // initialize time step variable
    int step = 0;

    // Save the macroscopic at t=0
    OutputVTK(step, this);
    OutputKeEns(step, this);

    // Simulation loop
    for (step = 1; step <= nstep; ++step)
    {
        #ifdef MULTICOMP
            #ifndef FD
            Collide_Species();  // collide species distribution function
            #else
            FD_species();       // finite difference for solving conservation of species
            #endif
        std::cout << "  Species Collision Done" << std::endl;
        #endif

        Collide();          // collision step
        std::cout << "  Mixture Collision Done" << std::endl;
        Streaming();        // streaming step & BC
        std::cout << "  Streaming Done" << std::endl;
        calculate_moment(); // calculate moment
        std::cout << "  Calculate Moment Done" << std::endl;

        if (step % tout == 0)
        {
            OutputVTK(step, this); // Save the macroscopic quantity
            OutputKeEns(step, this);
        }

    }
}


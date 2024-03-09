
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
                            species[a][i][j][k].rho = species[a][i][j][k].rho_n; 
                            Y[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho; 
                        }
                        gas->setMassFractions(&Y[0]);
                        // std::cout << 1.0/units.si_rho(mixture[i][j][k].rho) << std::endl;
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

void LBM::calculate_feq_geq(double f_tgt[], double g_tgt[], double rho_bb, double Y_bb[], double vel_tgt[], double temp_tgt)
{
    #ifndef MULTICOMP
    double p_tgt = rho_bb * gas_const * temp_tgt;
    double theta = gas_const*temp_tgt;
    double cv = gas_const / (gamma - 1.0);
    double cp = cv + gas_const;
    double internal_energy = cv * temp_tgt;
    double enthalpy = cp * temp_tgt; // H = Cp * T = (Cv + 1) * T
    #else
        int rank = omp_get_thread_num();
        auto gas = sols[rank]->thermo();
        std::vector <double> Y (gas->nSpecies());
        for(size_t a = 0; a < nSpecies; ++a)
            Y[gas->speciesIndex(speciesName[a])] = Y_bb[a];
        
        gas->setMassFractions(&Y[0]);
        gas->setState_TD(units.si_temp(temp_tgt), units.si_rho(rho_bb));

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
    }
}

// Initialize the disribution function from known macroscopic properties
void LBM::Init(){

    // Create Cantera's Solution object
    #ifdef MULTICOMP
    int nThreads = omp_get_max_threads();
    std::cout << "nthread : " << nThreads << std::endl;
    for(int i = 0; i < nThreads; ++i){
        if (diffModel == 0){
            auto sol = Cantera::newSolution("gri30.yaml", "gri30", "multicomponent");
            sols.push_back(sol);
        }
        else{ 
            auto sol = Cantera::newSolution("gri30.yaml", "gri30", "mixture-averaged");
            sols.push_back(sol);
        }
    }
    #endif
    

    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){    
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_I || mixture[i][j][k].type == TYPE_O){            
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

                    #ifdef MULTICOMP
                    // Species distribution function Initialization
                    for(size_t a = 0; a < nSpecies; ++a){
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
                if (mixture[i][j][k].type == TYPE_F){   
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

                    #ifdef MULTICOMP
                    double delYx[nSpecies], delYy[nSpecies], delYz[nSpecies];
                    for(size_t a = 0; a < nSpecies; ++a){
                        delYx[a] = ((species[a][i+1][j][k].rho / mixture[i+1][j][k].rho) - (species[a][i-1][j][k].rho / mixture[i-1][j][k].rho)) / (dx);
                        delYy[a] = ((species[a][i][j+1][k].rho / mixture[i][j+1][k].rho) - (species[a][i][j-1][k].rho / mixture[i][j-1][k].rho)) / (dy);
                        delYz[a] = ((species[a][i][j][k+1].rho / mixture[i][j][k+1].rho) - (species[a][i][j][k-1].rho / mixture[i][j][k-1].rho)) / (dz);
                    }
                    #endif

                    double delQdevx, delQdevy, delQdevz;
                    delQdevx = FD_limiterVanleer( mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1-3*gas_const*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1-3*gas_const*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u), dx) ;
                    delQdevy = FD_limiterVanleer( mixture[i][j-1][k].rho*mixture[i][j-1][k].u*(1-3*gas_const*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i][j+1][k].rho*mixture[i][j+1][k].u*(1-3*gas_const*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].u), dy) ;
                    delQdevz = FD_limiterVanleer( mixture[i][j][k-1].rho*mixture[i][j][k-1].u*(1-3*gas_const*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i][j][k+1].rho*mixture[i][j][k+1].u*(1-3*gas_const*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].u), dz) ;

                    double q_diff [3] = {0.0, 0.0, 0.0};
                    double q_corr [3] = {0.0, 0.0, 0.0};
                    #ifdef MULTICOMP
                    for (size_t a = 0; a < nSpecies; ++a){
                        int speciesIdx = gas->speciesIndex(speciesName[a]);
                        double mmass = gas->molecularWeight(speciesIdx);
                        q_diff[0] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].Vdiff_x);
                        q_diff[1] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].Vdiff_y);
                        q_diff[2] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].Vdiff_z);
                    
                        q_corr[0] += 0.5 * (omega1-2.0)/(omega1-omega) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * delYx[a];
                        q_corr[1] += 0.5 * (omega1-2.0)/(omega1-omega) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * delYy[a];
                        q_corr[2] += 0.5 * (omega1-2.0)/(omega1-omega) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * delYz[a];              
                    }
                    #endif

                    double str_heat_flux[3] ={  mixture[i][j][k].energy_flux[0] - velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) - velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) - velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) - 0.5*dt_sim*velocity[0]*delQdevx + q_diff[0] + q_corr[0],   //
                                                mixture[i][j][k].energy_flux[1] - velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) - velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) - velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) - 0.5*dt_sim*velocity[1]*delQdevy + q_diff[1] + q_corr[1],   //
                                                mixture[i][j][k].energy_flux[2] - velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) - velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) - velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) - 0.5*dt_sim*velocity[2]*delQdevz + q_diff[2] + q_corr[2]} ; //
                    
                    double corr[3] = {0.0, 0.0, 0.0};

                    // double corr[3] = {  dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevx,
                    //                     dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevy,
                    //                     dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevz};

                    for (int l = 0; l < npop; ++l){
                        // ------------- Mass and Momentum collision -----------------------------
                        double feq = calculate_feq(l, mixture[i][j][k].rho, velocity, theta, corr);
                        double geq = calculate_geq(l, mixture[i][j][k].rhoe, eq_heat_flux, eq_R_tensor, theta);
                        double gstr = calculate_geq(l, mixture[i][j][k].rhoe, str_heat_flux, eq_R_tensor, theta);

                        mixture[i][j][k].fpc[l] = mixture[i][j][k].f[l] + omega*(feq-mixture[i][j][k].f[l]);
                        mixture[i][j][k].gpc[l] = mixture[i][j][k].g[l] + omega1*(geq-mixture[i][j][k].g[l]) + (omega-omega1)*(gstr-mixture[i][j][k].g[l]);
                    }     
                }
            }
        }
    }
}

void LBM::Streaming(){
    // fill_FPC();

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
                        }
                        //---- Adiabatic Free-Slip Wall --------------------------
                        else if (mixture[i_nb][j_nb][k_nb].type==TYPE_FS){
                            int lp, ip, jp, kp;
                            dirSlip(l, i, j, k, lp, ip, jp, kp);
                            // std::cout << i << " | " << j << " | " << k << " | " << ip << " | " << jp << " | " << kp <<  std::endl;

                            mixture[i][j][k].f[l] = mixture[ip][jp][kp].fpc[lp];
                            mixture[i][j][k].g[l] = mixture[ip][jp][kp].gpc[lp];
                        }
                        //---- Adiabatic Wall --------------------------
                        else if (mixture[i_nb][j_nb][k_nb].type==TYPE_A){
                            mixture[i][j][k].f[l] = mixture[i][j][k].fpc[opposite[l]];
                            mixture[i][j][k].g[l] = mixture[i][j][k].gpc[opposite[l]];
                        }
                        //---- Inlet/Outlet Boundary Condition ---------------
                        else if (mixture[i_nb][j_nb][k_nb].type==TYPE_O){
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
                        }
                    }
                }
            }
        }
    }

    // // TMS Boundary Conditions ---------------------------------------------------------------------------------------------
    // // SOLID
    // #ifdef PARALLEL 
    // #pragma omp parallel for schedule(static, 1) 
    // #endif
    // for(int i=0; i<Nx; ++i){
    //     for(int j=0; j<Ny; ++j){
    //         for(int k = 0; k<Nz; ++k){
    //             if(mixture[i][j][k].type==TYPE_F){
    //                 int i_nb, j_nb, k_nb;
    //                 double f_tgt[npop];
    //                 double g_tgt[npop];
    //                 double f_loc[npop];
    //                 double g_loc[npop];
    //                 bool interface_nodes[npop] = {false}; 
    //                 int n_d = 0;

    //                 // check interface node                         
    //                 for (int l=0; l < npop; ++l){
    //                     i_nb = i-cx[l];
    //                     j_nb = j - cy[l];
    //                     k_nb = k - cz[l];

    //                     if(mixture[i_nb][j_nb][k_nb].type==TYPE_S && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
    //                         interface_nodes[l] = true;
    //                         n_d = n_d + 1;
    //                     }                            
    //                 }

    //                 // check if there is the node is the interface with boundary
    //                 if (n_d == 0)
    //                     continue;

    //                 // step 1: calculate f_tgt, g_tft, and fa_tgt
    //                 double rho_bb = 0.0;
    //                 double Y_bb[nSpecies] = {0.0};
    //                 for (int l=0; l < npop; ++l){
    //                     rho_bb += mixture[i][j][k].f[l];
    //                 }

    //                 double vel_tgt[3] = {0.0};
    //                 double T_tgt = 0.0;
    //                 for (int l=0; l < npop; ++l){
    //                     double q = 0.5;
    //                     if (interface_nodes[l] == true){
    //                         vel_tgt[0] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u)/(1.0+q);
    //                         vel_tgt[1] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].v+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v)/(1.0+q);
    //                         vel_tgt[2] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].w+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w)/(1.0+q);
    //                         T_tgt += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].temp+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp)/(1.0+q);
    //                         // std::cout << mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u << " | " << mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u << " | " << vel_tgt[0] << std::endl;
    //                     }
    //                 }

    //                 vel_tgt[0] = 1.0/n_d * vel_tgt[0];
    //                 vel_tgt[1] = 1.0/n_d * vel_tgt[1];
    //                 vel_tgt[2] = 1.0/n_d * vel_tgt[2];
    //                 T_tgt = 1.0/n_d * T_tgt;

    //                 // std::cout << vel_tgt[0] << " | " << T_tgt << std::endl;
    //                 calculate_feq_geq(f_tgt, g_tgt, rho_bb, Y_bb, vel_tgt, T_tgt);

    //                 // // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
    //                 double rho_loc = 0.0;
    //                 double rhou_loc = 0.0;
    //                 double rhov_loc = 0.0;
    //                 double rhow_loc = 0.0;
    //                 double rhoe_loc = 0.0;

    //                 for (int l=0; l < npop; ++l){
    //                     i_nb = i - cx[l];
    //                     j_nb = j - cy[l];
    //                     k_nb = k - cz[l];

    //                     if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
    //                         rho_loc += f_tgt[l];
    //                         rhou_loc += f_tgt[l]*cx[l];
    //                         rhov_loc += f_tgt[l]*cy[l];
    //                         rhow_loc += f_tgt[l]*cz[l];
    //                         rhoe_loc += g_tgt[l];
    //                     }
    //                     else{
    //                         rho_loc += mixture[i][j][k].f[l];
    //                         rhou_loc += mixture[i][j][k].f[l]*cx[l];
    //                         rhov_loc += mixture[i][j][k].f[l]*cy[l];
    //                         rhow_loc += mixture[i][j][k].f[l]*cz[l];
    //                         rhoe_loc += mixture[i][j][k].g[l];
    //                     }
    //                 }

    //                 double vel_loc[3];
    //                 vel_loc[0] = rhou_loc / rho_loc;
    //                 vel_loc[1] = rhov_loc / rho_loc;
    //                 vel_loc[2] = rhow_loc / rho_loc;

    //                 double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                    
    //                 #ifndef MULTICOMP
    //                     double cv = gas_const / (gamma - 1.0);
    //                     double T_loc = internalEnergy / cv;                        
                        
    //                 #else
    //                     int rank = omp_get_thread_num();
    //                     auto gas = sols[rank]->thermo();
    //                     std::vector <double> Y_loc (gas->nSpecies());

    //                     for(size_t a = 0; a < nSpecies; ++a){
    //                         Y_loc[gas->speciesIndex(speciesName[a])] = rhoa_loc[a] / rho_loc;
    //                         std::cout << rhoa_loc[a] / rho_loc << std::endl;
    //                     }
    //                     gas->setMassFractions(&Y_loc[0]);
    //                     gas->setState_UV(units.si_energy_mass(internalEnergy), 1.0 / units.si_rho(rho_loc));

    //                     double T_loc = gas->temperature();
    //                 #endif

    //                 calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);

    //                 for (int l=0; l < npop; ++l){
    //                     i_nb = i - cx[l];
    //                     j_nb = j - cy[l];
    //                     k_nb = k - cz[l];

    //                     if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
    //                         mixture[i][j][k].f[l] = 2*f_tgt[l] - f_loc[l];
    //                         mixture[i][j][k].g[l] = 2*g_tgt[l] - g_loc[l];
    //                     }
    //                     else{
    //                         mixture[i][j][k].f[l] = f_tgt[l] + mixture[i][j][k].f[l] - f_loc[l];
    //                         mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
    //                     }
    //                 }                   

    //             } 
    //         }
    //     }
    // }

    // // INFLOW
    // #ifdef PARALLEL 
    // #pragma omp parallel for schedule(static, 1) 
    // #endif
    // for(int i=0; i<Nx; ++i){
    //     for(int j=0; j<Ny; ++j){
    //         for(int k = 0; k<Nz; ++k){
    //             if(mixture[i][j][k].type==TYPE_F){
    //                 int i_nb, j_nb, k_nb;
    //                 double f_in[npop];
    //                 double g_in[npop];
    //                 double f_tgt[npop];
    //                 double g_tgt[npop];
    //                 double f_loc[npop];
    //                 double g_loc[npop];
    //                 bool interface_nodes[npop] = {false}; 
    //                 int n_d = 0;

    //                 // check interface node                         
    //                 for (int l=0; l < npop; ++l){
    //                     i_nb = i - cx[l];
    //                     j_nb = j - cy[l];
    //                     k_nb = k - cz[l];

    //                     if(mixture[i_nb][j_nb][k_nb].type==TYPE_I && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
    //                         interface_nodes[l] = true;
    //                         n_d = n_d + 1;
    //                         break;
    //                     }                            
    //                 }

    //                 // check if there is the node is the interface with boundary
    //                 if (n_d == 0)
    //                     continue;

    //                 double vel_in[3] = {0.0};
    //                 double T_in = 0.0;
    //                 double rho_in = 0.0;
    //                 for (int l=0; l < npop; ++l){
    //                     if (interface_nodes[l] == true){
    //                         vel_in[0] = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u;
    //                         vel_in[1] = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v;
    //                         vel_in[2] = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w;
    //                         T_in = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp;
    //                         rho_in = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].rho;
    //                         break;
    //                     }
    //                 }

    //                 // std::cout << i << " | " << vel_in[0] << " | " << T_in << std::endl;
    //                 calculate_feq_geq(f_in, g_in, rho_in, vel_in, T_in);

    //                 // // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
    //                 double rho_loc = 0.0;
    //                 double rhou_loc = 0.0;
    //                 double rhov_loc = 0.0;
    //                 double rhow_loc = 0.0;
    //                 double rhoe_loc = 0.0;

    //                 for (int l=0; l < npop; ++l){
    //                     i_nb = i - cx[l];
    //                     j_nb = j - cy[l];
    //                     k_nb = k - cz[l];

    //                     if (mixture[i_nb][j_nb][k_nb].type==TYPE_I){
    //                         rho_loc += f_in[l];
    //                         rhou_loc += f_in[l]*cx[l];
    //                         rhov_loc += f_in[l]*cy[l];
    //                         rhow_loc += f_in[l]*cz[l];
    //                         rhoe_loc += g_in[l];
    //                     }
    //                     else{
    //                         rho_loc += mixture[i][j][k].f[l];
    //                         rhou_loc += mixture[i][j][k].f[l]*cx[l];
    //                         rhov_loc += mixture[i][j][k].f[l]*cy[l];
    //                         rhow_loc += mixture[i][j][k].f[l]*cz[l];
    //                         rhoe_loc += mixture[i][j][k].g[l];
    //                     }
    //                 }

    //                 double vel_loc[3];
    //                 vel_loc[0] = rhou_loc / rho_loc;
    //                 vel_loc[1] = rhov_loc / rho_loc;
    //                 vel_loc[2] = rhow_loc / rho_loc;

    //                 double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                    
    //                 #ifndef MULTICOMP
    //                     double cv = gas_const / (gamma - 1.0);
    //                     double T_loc = internalEnergy / cv;                        
                        
    //                 #else
    //                     int rank = omp_get_thread_num();
    //                     auto gas = sols[rank]->thermo();
    //                     std::vector <double> Y_loc (gas->nSpecies());

    //                     for(size_t a = 0; a < nSpecies; ++a){
    //                         Y_loc[gas->speciesIndex(speciesName[a])] = rhoa_loc[a] / rho_loc;
    //                         std::cout << rhoa_loc[a] / rho_loc << std::endl;
    //                     }
    //                     gas->setMassFractions(&Y_loc[0]);
    //                     gas->setState_UV(units.si_energy_mass(internalEnergy), 1.0 / units.si_rho(rho_loc));

    //                     double T_loc = gas->temperature();
    //                 #endif

    //                 calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);

    //                 calculate_feq_geq(f_tgt, g_tgt, rho_loc, vel_in, T_in);


    //                 for (int l=0; l < npop; ++l){
    //                     i_nb = i - cx[l];
    //                     j_nb = j - cy[l];
    //                     k_nb = k - cz[l];

    //                     if (mixture[i_nb][j_nb][k_nb].type==TYPE_I){
    //                         mixture[i][j][k].f[l] = f_tgt[l] + f_in[l] - f_loc[l];
    //                         mixture[i][j][k].g[l] = g_tgt[l] + g_in[l] - g_loc[l];
    //                     }
    //                     else{
    //                         mixture[i][j][k].f[l] = f_tgt[l] + mixture[i][j][k].f[l] - f_loc[l];
    //                         mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
    //                     }
    //                 }                   

    //             } 
    //         }
    //     }
    // }

    // // OUTFLOW
    // #ifdef PARALLEL 
    // #pragma omp parallel for schedule(static, 1) 
    // #endif
    // for(int i=0; i<Nx; ++i){
    //     for(int j=0; j<Ny; ++j){
    //         for(int k = 0; k<Nz; ++k){
    //             if(mixture[i][j][k].type==TYPE_F){
    //                 int i_nb, j_nb, k_nb;
    //                 double f_out[npop];
    //                 double g_out[npop];
    //                 double f_loc[npop];
    //                 double g_loc[npop];
    //                 bool interface_nodes[npop] = {false}; 
    //                 int n_d = 0;

    //                 // check interface node                         
    //                 for (int l=0; l < npop; ++l){
    //                     i_nb = i - cx[l];
    //                     j_nb = j - cy[l];
    //                     k_nb = k - cz[l];

    //                     if(mixture[i_nb][j_nb][k_nb].type==TYPE_O && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
    //                         interface_nodes[l] = true;
    //                         n_d = n_d + 1;
    //                         break;
    //                     }                            
    //                 }

    //                 // check if there is the node is the interface with boundary
    //                 if (n_d == 0)
    //                     continue;

    //                 double vel_out[3] = {0.0};
    //                 double T_out = 0.0;
    //                 double rho_out = 0.0;
    //                 for (int l=0; l < npop; ++l){
    //                     if (interface_nodes[l] == true){
    //                         vel_out[0] = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u;
    //                         vel_out[1] = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v;
    //                         vel_out[2] = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w;
    //                         T_out = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp;
    //                         rho_out = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].rho;
    //                         break;
    //                     }
    //                 }

    //                 // std::cout << i << " | " << vel_out[0] << " | " << T_out << std::endl;
    //                 calculate_feq_geq(f_out, g_out, rho_out, vel_out, T_out);

    //                 // // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
    //                 double rho_loc = 0.0;
    //                 double rhou_loc = 0.0;
    //                 double rhov_loc = 0.0;
    //                 double rhow_loc = 0.0;
    //                 double rhoe_loc = 0.0;

    //                 for (int l=0; l < npop; ++l){
    //                     i_nb = i - cx[l];
    //                     j_nb = j - cy[l];
    //                     k_nb = k - cz[l];

    //                     if (mixture[i_nb][j_nb][k_nb].type==TYPE_O){
    //                         rho_loc += f_out[l];
    //                         rhou_loc += f_out[l]*cx[l];
    //                         rhov_loc += f_out[l]*cy[l];
    //                         rhow_loc += f_out[l]*cz[l];
    //                         rhoe_loc += g_out[l];
    //                     }
    //                     else{
    //                         rho_loc += mixture[i][j][k].f[l];
    //                         rhou_loc += mixture[i][j][k].f[l]*cx[l];
    //                         rhov_loc += mixture[i][j][k].f[l]*cy[l];
    //                         rhow_loc += mixture[i][j][k].f[l]*cz[l];
    //                         rhoe_loc += mixture[i][j][k].g[l];
    //                     }
    //                 }

    //                 double vel_loc[3];
    //                 vel_loc[0] = rhou_loc / rho_loc;
    //                 vel_loc[1] = rhov_loc / rho_loc;
    //                 vel_loc[2] = rhow_loc / rho_loc;

    //                 double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                    
    //                 #ifndef MULTICOMP
    //                     double cv = gas_const / (gamma - 1.0);
    //                     double T_loc = internalEnergy / cv;                        
                        
    //                 #else
    //                     int rank = omp_get_thread_num();
    //                     auto gas = sols[rank]->thermo();
    //                     std::vector <double> Y_loc (gas->nSpecies());

    //                     for(size_t a = 0; a < nSpecies; ++a){
    //                         Y_loc[gas->speciesIndex(speciesName[a])] = rhoa_loc[a] / rho_loc;
    //                         std::cout << rhoa_loc[a] / rho_loc << std::endl;
    //                     }
    //                     gas->setMassFractions(&Y_loc[0]);
    //                     gas->setState_UV(units.si_energy_mass(internalEnergy), 1.0 / units.si_rho(rho_loc));

    //                     double T_loc = gas->temperature();
    //                 #endif

    //                 calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);

    //                 for (int l=0; l < npop; ++l){
    //                     i_nb = i - cx[l];
    //                     j_nb = j - cy[l];
    //                     k_nb = k - cz[l];

    //                     if (mixture[i_nb][j_nb][k_nb].type==TYPE_O){
    //                         mixture[i][j][k].f[l] = 2.0*f_out[l] - f_loc[l];
    //                         mixture[i][j][k].g[l] = 2.0*g_out[l] - g_loc[l];
    //                     }
    //                     else{
    //                         mixture[i][j][k].f[l] = f_out[l] + mixture[i][j][k].f[l] - f_loc[l];
    //                         mixture[i][j][k].g[l] = g_out[l] + mixture[i][j][k].g[l] - g_loc[l];
    //                     }
    //                 }                   

    //             } 
    //         }
    //     }
    // }


}



                   


#ifdef MULTICOMP
void LBM::FD_species()
{   
    #ifdef PARALLEL 
    #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){    
                if (mixture[i][j][k].type == TYPE_F){
                    int rank = omp_get_thread_num();
                    auto gas = sols[rank]->thermo();               
                    std::vector<double> Y (gas->nSpecies());
                    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));
                    auto trans = sols[rank]->transport();
                    
                    std::vector<double> delT(ndim);
                    delT[0] = FD_limiterVanleer(mixture[i-1][j][k].temp, mixture[i][j][k].temp, mixture[i+1][j][k].temp, dx);
                    delT[1] = FD_limiterVanleer(mixture[i][j-1][k].temp, mixture[i][j][k].temp, mixture[i][j+1][k].temp, dy);
                    delT[2] = FD_limiterVanleer(mixture[i][j][k-1].temp, mixture[i][j][k].temp, mixture[i][j][k+1].temp, dz);
                    
                    std::vector<double> delX(ndim*nSpecies);
                    for (size_t a = 0; a < nSpecies; ++a)
                        delX[a]              = FD_limiterVanleer(species[a][i-1][j][k].X, species[a][i][j][k].X, species[a][i+1][j][k].X, dx);
                    for (size_t a = 0; a < nSpecies; ++a)
                        delX[nSpecies + a]   = FD_limiterVanleer(species[a][i][j-1][k].X, species[a][i][j][k].X, species[a][i][j+1][k].X, dy);
                    for (size_t a = 0; a < nSpecies; ++a)
                        delX[2*nSpecies + a] = FD_limiterVanleer(species[a][i][j][k-1].X, species[a][i][j][k].X, species[a][i][j][k+1].X, dz);
                    
                    std::vector<double> delP(ndim*nSpecies);
                    for (size_t a = 0; a < nSpecies; ++a)
                        delP[a]                 = FD_limiterVanleer(mixture[i-1][j][k].p, mixture[i][j][k].p, mixture[i+1][j][k].p, dx);
                    for (size_t a = 0; a < nSpecies; ++a)
                        delP[nSpecies + a]      = FD_limiterVanleer(mixture[i][j-1][k].p, mixture[i][j][k].p, mixture[i][j+1][k].p, dy);
                    for (size_t a = 0; a < nSpecies; ++a)
                        delP[2*nSpecies + a]    = FD_limiterVanleer(mixture[i][j][k-1].p, mixture[i][j][k].p, mixture[i][j][k+1].p, dz);

                    std::vector<double> diff_force(ndim*nSpecies);
                    for (size_t a = 0; a < nSpecies; ++a)
                        diff_force[a]               = delX[a] ;//+ (species[a][i][j][k].X - species[a][i][j][k].rho/mixture[i][j][k].rho) * delP[a]/mixture[i][j][k].p;
                    for (size_t a = 0; a < nSpecies; ++a)
                        diff_force[nSpecies + a]    = delX[nSpecies+a] ;//+ (species[a][i][j][k].X - species[a][i][j][k].rho/mixture[i][j][k].rho) * delP[nSpecies+a]/mixture[i][j][k].p;
                    for (size_t a = 0; a < nSpecies; ++a)
                        diff_force[2*nSpecies + a]  = delX[2*nSpecies+a] ;//+ (species[a][i][j][k].X - species[a][i][j][k].rho/mixture[i][j][k].rho) * delP[2*nSpecies+a]/mixture[i][j][k].p;

                    std::vector<double> thermal_diff_coeff(gas->nSpecies());
                    trans->getThermalDiffCoeffs(&thermal_diff_coeff[0]);

                    if (diffModel == 0){
                        std::vector<double> multi_diff_coeff(gas->nSpecies()*gas->nSpecies());
                        trans->getMultiDiffCoeffs(gas->nSpecies(), &multi_diff_coeff[0]);

                        double sum[ndim*nSpecies] = {0.0};
                        for (size_t a = 0; a < nSpecies; ++a){
                            for (size_t b = 0; b < nSpecies; ++b)
                                sum[a]              += gas->molecularWeight(gas->speciesIndex(speciesName[b]))*units.nu(multi_diff_coeff[gas->nSpecies()*gas->speciesIndex(speciesName[a])+gas->speciesIndex(speciesName[b])])*diff_force[b];
                        }
                        for (size_t a = 0; a < nSpecies; ++a){
                            for (size_t b = 0; b < nSpecies; ++b)
                                sum[nSpecies + a]   += gas->molecularWeight(gas->speciesIndex(speciesName[b]))*units.nu(multi_diff_coeff[gas->nSpecies()*gas->speciesIndex(speciesName[a])+gas->speciesIndex(speciesName[b])])*diff_force[nSpecies + b];
                        }
                        for (size_t a = 0; a < nSpecies; ++a){
                            for (size_t b = 0; b < nSpecies; ++b)
                                sum[2*nSpecies + a] += gas->molecularWeight(gas->speciesIndex(speciesName[b]))*units.nu(multi_diff_coeff[gas->nSpecies()*gas->speciesIndex(speciesName[a])+gas->speciesIndex(speciesName[b])])*diff_force[2*nSpecies + b];
                        }

                        for (size_t a = 0; a < nSpecies; ++a){
                            if (species[a][i][j][k].X >= 1e-20)
                                species[a][i][j][k].Vdiff_x = 1.0/species[a][i][j][k].X/gas->meanMolecularWeight()*sum[a]              ;//- units.mu(thermal_diff_coeff[gas->speciesIndex(speciesName[a])])/species[a][i][j][k].rho/mixture[i][j][k].temp * delT[0];
                            else
                                species[a][i][j][k].Vdiff_x = 0.0;
                        }
                        for (size_t a = 0; a < nSpecies; ++a){
                            if (species[a][i][j][k].X >= 1e-20)
                                species[a][i][j][k].Vdiff_y = 1.0/species[a][i][j][k].X/gas->meanMolecularWeight()*sum[nSpecies + a]   ;//- units.mu(thermal_diff_coeff[gas->speciesIndex(speciesName[a])])/species[a][i][j][k].rho/mixture[i][j][k].temp * delT[1];
                            else
                                species[a][i][j][k].Vdiff_y = 0.0;
                        }
                        for (size_t a = 0; a < nSpecies; ++a){
                            if (species[a][i][j][k].X >= 1e-20)
                                species[a][i][j][k].Vdiff_z = 1.0/species[a][i][j][k].X/gas->meanMolecularWeight()*sum[2*nSpecies + a] ;//- units.mu(thermal_diff_coeff[gas->speciesIndex(speciesName[a])])/species[a][i][j][k].rho/mixture[i][j][k].temp * delT[2];
                            else
                                species[a][i][j][k].Vdiff_z = 0.0;
                        }
                    }
                    else{
                        std::vector<double> mix_diff_coeff(gas->nSpecies());
                        trans->getMixDiffCoeffsMole(&mix_diff_coeff[0]);
                        for (size_t a = 0; a < nSpecies; ++a)
                            if (species[a][i][j][k].X >= 1e-20)
                                species[a][i][j][k].Vdiff_x = -1.0/species[a][i][j][k].X*units.nu(mix_diff_coeff[gas->speciesIndex(speciesName[a])])*diff_force[a] - 1.0/mixture[i][j][k].temp*units.mu(thermal_diff_coeff[gas->speciesIndex(speciesName[a])])/(species[a][i][j][k].rho)*delT[0];
                            else
                                species[a][i][j][k].Vdiff_x = 0.0;
                        for (size_t a = 0; a < nSpecies; ++a)
                            if (species[a][i][j][k].X >= 1e-20)
                                species[a][i][j][k].Vdiff_y = -1.0/species[a][i][j][k].X*units.nu(mix_diff_coeff[gas->speciesIndex(speciesName[a])])*diff_force[nSpecies+a] - 1.0/mixture[i][j][k].temp*units.mu(thermal_diff_coeff[gas->speciesIndex(speciesName[a])])/(species[a][i][j][k].rho)*delT[1];
                            else
                                species[a][i][j][k].Vdiff_y = 0.0;
                        for (size_t a = 0; a < nSpecies; ++a)
                            if (species[a][i][j][k].X >= 1e-20)
                                species[a][i][j][k].Vdiff_z = -1.0/species[a][i][j][k].X*units.nu(mix_diff_coeff[gas->speciesIndex(speciesName[a])])*diff_force[2*nSpecies+a] - 1.0/mixture[i][j][k].temp*units.mu(thermal_diff_coeff[gas->speciesIndex(speciesName[a])])/(species[a][i][j][k].rho)*delT[2];
                            else
                                species[a][i][j][k].Vdiff_z = 0.0;
                    }

                    double Vcorr_x = 0.0;
                    double Vcorr_y = 0.0;
                    double Vcorr_z = 0.0;
                    for (size_t a = 0; a < nSpecies; ++a){
                        Vcorr_x += species[a][i][j][k].rho / mixture[i][j][k].rho * species[a][i][j][k].Vdiff_x;
                        Vcorr_y += species[a][i][j][k].rho / mixture[i][j][k].rho * species[a][i][j][k].Vdiff_y;
                        Vcorr_z += species[a][i][j][k].rho / mixture[i][j][k].rho * species[a][i][j][k].Vdiff_z;
                    }

                    for (size_t a = 0; a < nSpecies; ++a){
                        species[a][i][j][k].Vdiff_x -= Vcorr_x;
                        species[a][i][j][k].Vdiff_y -= Vcorr_y;
                        species[a][i][j][k].Vdiff_z -= Vcorr_z;
                    }
                    // std::cout << i << " | " << species[0][i][j][k].Vdiff_x << " | " << species[1][i][j][k].Vdiff_x << std::endl;
                }
            }
        }
    }

    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){    
                if (mixture[i][j][k].type == TYPE_F){    
                    for (size_t a = 0; a < nSpecies; ++a){
                        // double conv_x = FD_limiterVanleer(species[a][i-1][j][k].rho*(mixture[i-1][j][k].u + species[a][i-1][j][k].Vdiff_x) , species[a][i][j][k].rho*(mixture[i][j][k].u + species[a][i][j][k].Vdiff_x), species[a][i+1][j][k].rho*(mixture[i+1][j][k].u + species[a][i+1][j][k].Vdiff_x), dx );
                        // double conv_y = FD_limiterVanleer(species[a][i][j-1][k].rho*(mixture[i][j-1][k].v + species[a][i][j-1][k].Vdiff_y) , species[a][i][j][k].rho*(mixture[i][j][k].v + species[a][i][j][k].Vdiff_y), species[a][i][j+1][k].rho*(mixture[i][j+1][k].v + species[a][i][j+1][k].Vdiff_y), dy );
                        // double conv_z = FD_limiterVanleer(species[a][i][j][k-1].rho*(mixture[i][j][k-1].w + species[a][i][j][k-1].Vdiff_z) , species[a][i][j][k].rho*(mixture[i][j][k].w + species[a][i][j][k].Vdiff_z), species[a][i][j][k+1].rho*(mixture[i][j][k+1].w + species[a][i][j][k+1].Vdiff_z), dz );

                        double conv_x = (species[a][i+1][j][k].rho*(mixture[i+1][j][k].u + species[a][i+1][j][k].Vdiff_x) - species[a][i-1][j][k].rho*(mixture[i-1][j][k].u + species[a][i-1][j][k].Vdiff_x)) / (2*dx);
                        double conv_y = (species[a][i][j+1][k].rho*(mixture[i][j+1][k].v + species[a][i][j+1][k].Vdiff_y) - species[a][i][j-1][k].rho*(mixture[i][j-1][k].v + species[a][i][j-1][k].Vdiff_y)) / (2*dy);
                        double conv_z = (species[a][i][j][k+1].rho*(mixture[i][j][k+1].w + species[a][i][j][k+1].Vdiff_z) - species[a][i][j][k-1].rho*(mixture[i][j][k-1].w + species[a][i][j][k-1].Vdiff_z)) / (2*dz);

                        // std::cout << i << " | " << conv_x << " | " << species[a][i-1][j][k].rho*(mixture[i-1][j][k].u + species[a][i-1][j][k].Vdiff_x) << " | " << species[a][i][j][k].rho*(mixture[i][j][k].u + species[a][i][j][k].Vdiff_x) << " | " << species[a][i][j+1][k].rho*(mixture[i+1][j][k].u + species[a][i+1][j][k].Vdiff_x) << std::endl;
                        species[a][i][j][k].rho_n = species[a][i][j][k].rho - dt_sim * ( conv_x + conv_y + conv_z );           
                    }

                }
            }
        }
    }
}
#endif

void LBM::dirSlip(int l, int i, int j, int k, int &lp, int &ip, int &jp, int &kp)
{
    int i_nb, j_nb, k_nb;
    i_nb = i - cx[l];
    j_nb = j - cy[l];
    k_nb = k - cz[l];

    bool i_chk = false;
    bool j_chk = false;
    bool k_chk = false;
    if (mixture[i_nb][j][k].type==TYPE_FS)
        i_chk = true;
    if (mixture[i][j_nb][k].type==TYPE_FS)
        j_chk = true;
    if (mixture[i][j][k_nb].type==TYPE_FS)
        k_chk = true;

    double cxp = cx[l];
    double cyp = cy[l];
    double czp = cz[l];
    ip = i_nb;
    jp = j_nb;
    kp = k_nb;
    if (i_chk == true){
        cxp = -1.0*cx[l];
        ip = i;
    }
    if (j_chk == true){
        cyp = -1.0*cy[l];
        jp = j;
    }
    if (k_chk == true){
        czp = -1.0*cz[l];
        kp = k;
    } 
    
    for(size_t a = 0; a < npop; ++a){
        if (cxp == cx[a] && cyp == cy[a] && czp == cz[a]){
            lp = a;
            break;
        }
    }

}

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
                    if(mixture[i_nb][j][k].type == TYPE_A || mixture[i_nb][j][k].type == TYPE_FS){
                        mixture[i_nb][j][k].rho = mixture[i][j][k].rho;
                        mixture[i_nb][j][k].u   = mixture[i][j][k].u;
                        mixture[i_nb][j][k].v   = mixture[i][j][k].v;
                        mixture[i_nb][j][k].w   = mixture[i][j][k].w;
                        mixture[i_nb][j][k].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_nb][j][k].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i_pb][j][k].type == TYPE_A || mixture[i_pb][j][k].type == TYPE_FS){
                        mixture[i_pb][j][k].rho = mixture[i][j][k].rho;
                        mixture[i_pb][j][k].u   = mixture[i][j][k].u;
                        mixture[i_pb][j][k].v   = mixture[i][j][k].v;
                        mixture[i_pb][j][k].w   = mixture[i][j][k].w;
                        mixture[i_pb][j][k].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_pb][j][k].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j_nb][k].type == TYPE_A || mixture[i][j_nb][k].type == TYPE_FS){
                        mixture[i][j_nb][k].rho = mixture[i][j][k].rho;
                        mixture[i][j_nb][k].u   = mixture[i][j][k].u;
                        mixture[i][j_nb][k].v   = mixture[i][j][k].v;
                        mixture[i][j_nb][k].w   = mixture[i][j][k].w;
                        mixture[i][j_nb][k].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_nb][k].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j_pb][k].type == TYPE_A || mixture[i][j_pb][k].type == TYPE_FS){
                        mixture[i][j_pb][k].rho = mixture[i][j][k].rho;
                        mixture[i][j_pb][k].u   = mixture[i][j][k].u;
                        mixture[i][j_pb][k].v   = mixture[i][j][k].v;
                        mixture[i][j_pb][k].w   = mixture[i][j][k].w;
                        mixture[i][j_pb][k].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_pb][k].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j][k_nb].type == TYPE_A || mixture[i][j][k_nb].type == TYPE_FS){
                        mixture[i][j][k_nb].rho = mixture[i][j][k].rho;
                        mixture[i][j][k_nb].u   = mixture[i][j][k].u;
                        mixture[i][j][k_nb].v   = mixture[i][j][k].v;
                        mixture[i][j][k_nb].w   = mixture[i][j][k].w;
                        mixture[i][j][k_nb].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_nb].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j][k_pb].type == TYPE_A || mixture[i][j][k_pb].type == TYPE_FS){
                        mixture[i][j][k_pb].rho = mixture[i][j][k].rho;
                        mixture[i][j][k_pb].u   = mixture[i][j][k].u;
                        mixture[i][j][k_pb].v   = mixture[i][j][k].v;
                        mixture[i][j][k_pb].w   = mixture[i][j][k].w;
                        mixture[i][j][k_pb].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_pb].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    
                    // Periodic Boundary Condition 
                    if(mixture[i_nb][j][k].type == TYPE_P){
                        mixture[i_nb][j][k].rho = mixture[Nx-2][j][k].rho;
                        mixture[i_nb][j][k].u   = mixture[Nx-2][j][k].u;
                        mixture[i_nb][j][k].v   = mixture[Nx-2][j][k].v;
                        mixture[i_nb][j][k].w   = mixture[Nx-2][j][k].w;
                        mixture[i_nb][j][k].temp= mixture[Nx-2][j][k].temp;
                        for(int l = 0; l < npop; ++l){
                            mixture[i_nb][j][k].f[l] = mixture[Nx-2][j][k].f[l];
                            mixture[i_nb][j][k].g[l] = mixture[Nx-2][j][k].g[l];
                        }
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_nb][j][k].rho = species[a][Nx-2][j][k].rho;
                        #endif
                    }
                    if(mixture[i_pb][j][k].type == TYPE_P){
                        mixture[i_pb][j][k].rho = mixture[1][j][k].rho;
                        mixture[i_pb][j][k].u   = mixture[1][j][k].u;
                        mixture[i_pb][j][k].v   = mixture[1][j][k].v;
                        mixture[i_pb][j][k].w   = mixture[1][j][k].w;
                        mixture[i_pb][j][k].temp= mixture[1][j][k].temp;
                        for(int l = 0; l < npop; ++l){
                            mixture[i_pb][j][k].f[l] = mixture[1][j][k].f[l];
                            mixture[i_pb][j][k].g[l] = mixture[1][j][k].g[l];
                        }
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_pb][j][k].rho = species[a][1][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j_nb][k].type == TYPE_P){
                        mixture[i][j_nb][k].rho = mixture[i][Ny-2][k].rho;
                        mixture[i][j_nb][k].u   = mixture[i][Ny-2][k].u;
                        mixture[i][j_nb][k].v   = mixture[i][Ny-2][k].v;
                        mixture[i][j_nb][k].w   = mixture[i][Ny-2][k].w;
                        mixture[i][j_nb][k].temp= mixture[i][Ny-2][k].temp;
                        for(int l = 0; l < npop; ++l){
                            mixture[i][j_nb][k].f[l] = mixture[i][Ny-2][k].f[l];
                            mixture[i][j_nb][k].g[l] = mixture[i][Ny-2][k].g[l];
                        }
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_nb][k].rho = species[a][i][Ny-2][k].rho;
                        #endif
                    }
                    if(mixture[i][j_pb][k].type == TYPE_P){
                        mixture[i][j_pb][k].rho = mixture[i][1][k].rho;
                        mixture[i][j_pb][k].u   = mixture[i][1][k].u;
                        mixture[i][j_pb][k].v   = mixture[i][1][k].v;
                        mixture[i][j_pb][k].w   = mixture[i][1][k].w;
                        mixture[i][j_pb][k].temp= mixture[i][1][k].temp;
                        for(int l = 0; l < npop; ++l){
                            mixture[i][j_pb][k].f[l] = mixture[i][1][k].f[l];
                            mixture[i][j_pb][k].g[l] = mixture[i][1][k].g[l];
                        }
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_pb][k].rho = species[a][i][1][k].rho;
                        #endif
                    }
                    if(mixture[i][j][k_nb].type == TYPE_P){
                        mixture[i][j][k_nb].rho = mixture[i][j][Nz-2].rho;
                        mixture[i][j][k_nb].u   = mixture[i][j][Nz-2].u;
                        mixture[i][j][k_nb].v   = mixture[i][j][Nz-2].v;
                        mixture[i][j][k_nb].w   = mixture[i][j][Nz-2].w;
                        mixture[i][j][k_nb].temp= mixture[i][j][Nz-2].temp;
                        for(int l = 0; l < npop; ++l){
                            mixture[i][j][k_nb].f[l] = mixture[i][j][Nz-2].f[l];
                            mixture[i][j][k_nb].g[l] = mixture[i][j][Nz-2].g[l];
                        }
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_nb].rho = species[a][i][j][Nz-2].rho;
                        #endif
                    }
                    if(mixture[i][j][k_pb].type == TYPE_P){
                        mixture[i][j][k_pb].rho = mixture[i][j][1].rho;
                        mixture[i][j][k_pb].u   = mixture[i][j][1].u;
                        mixture[i][j][k_pb].v   = mixture[i][j][1].v;
                        mixture[i][j][k_pb].w   = mixture[i][j][1].w;
                        mixture[i][j][k_pb].temp= mixture[i][j][1].temp;
                        for(int l = 0; l < npop; ++l){
                            mixture[i][j][k_pb].f[l] = mixture[i][j][1].f[l];
                            mixture[i][j][k_pb].g[l] = mixture[i][j][1].g[l];
                        }
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_pb].rho = species[a][i][j][1].rho;
                        #endif
                    }


                }
            }
        }
    }
}

void LBM::fill_FPC()
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

                                        
                    // Periodic Boundary Condition 
                    if(mixture[i_nb][j][k].type == TYPE_P){
                        for(int l = 0; l < npop; ++l){
                            mixture[i_nb][j][k].fpc[l] = mixture[Nx-2][j][k].fpc[l];
                            mixture[i_nb][j][k].gpc[l] = mixture[Nx-2][j][k].gpc[l];
                        }
                    }
                    if(mixture[i_pb][j][k].type == TYPE_P){
                        for(int l = 0; l < npop; ++l){
                            mixture[i_pb][j][k].fpc[l] = mixture[1][j][k].fpc[l];
                            mixture[i_pb][j][k].gpc[l] = mixture[1][j][k].gpc[l];
                        }
                    }
                    if(mixture[i][j_nb][k].type == TYPE_P){
                        for(int l = 0; l < npop; ++l){
                            mixture[i][j_nb][k].fpc[l] = mixture[i][Ny-2][k].fpc[l];
                            mixture[i][j_nb][k].gpc[l] = mixture[i][Ny-2][k].gpc[l];
                        }
                    }
                    if(mixture[i][j_pb][k].type == TYPE_P){
                        for(int l = 0; l < npop; ++l){
                            mixture[i][j_pb][k].fpc[l] = mixture[i][1][k].fpc[l];
                            mixture[i][j_pb][k].gpc[l] = mixture[i][1][k].gpc[l];
                        }
                    }
                    if(mixture[i][j][k_nb].type == TYPE_P){
                        for(int l = 0; l < npop; ++l){
                            mixture[i][j][k_nb].fpc[l] = mixture[i][j][Nz-2].fpc[l];
                            mixture[i][j][k_nb].gpc[l] = mixture[i][j][Nz-2].gpc[l];
                        }
                    }
                    if(mixture[i][j][k_pb].type == TYPE_P){
                        for(int l = 0; l < npop; ++l){
                            mixture[i][j][k_pb].fpc[l] = mixture[i][j][1].fpc[l];
                            mixture[i][j][k_pb].gpc[l] = mixture[i][j][1].gpc[l];
                        }
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
        FD_species();       // finite difference for solving conservation of species
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


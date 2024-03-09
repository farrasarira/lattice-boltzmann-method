
#include "cantera.hpp"
#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include <omp.h>
#include <numeric>
#include <vector>

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

void LBM::calculate_moment()
{
    std::vector<std::vector<std::vector<double>>> dQdevx(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevy(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    std::vector<std::vector<std::vector<double>>> dQdevz(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));

    #ifdef PARALLEL 
    #pragma omp parallel for schedule(static, 1) 
    #endif
    for (int i = 0; i < Nx; ++i){
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

                    double cv = gas_const / (gamma - 1.0);
                    mixture[i][j][k].temp = internalEnergy / cv;
                    mixture[i][j][k].p = mixture[i][j][k].rho*gas_const*mixture[i][j][k].temp;
            
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

// Initialize the disribution function from known macroscopic properties
void LBM::Init(){

    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){    
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_I || mixture[i][j][k].type == TYPE_O){               
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
                    double cv = gas_const / (gamma - 1.0);
                    double cp = cv + gas_const;
                    double theta = gas_const*mixture[i][j][k].temp;
                    double mu = nu*mixture[i][j][k].rho;
                    double conduc_coeff = mu*cp/prtl;
                    // std::cout << "alpha : " << conduc_coeff/(cp * mixture[i][j][k].rho) << std::endl;
                    double enthalpy = cp * mixture[i][j][k].temp; // H = Cp * T = (Cv + 1) * T
                
                    

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

                    double delQdevx, delQdevy, delQdevz;
                    delQdevx = FD_limiterVanleer( mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1-3*gas_const*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1-3*gas_const*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u), dx) ;
                    delQdevy = FD_limiterVanleer( mixture[i][j-1][k].rho*mixture[i][j-1][k].u*(1-3*gas_const*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i][j+1][k].rho*mixture[i][j+1][k].u*(1-3*gas_const*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].u), dy) ;
                    delQdevz = FD_limiterVanleer( mixture[i][j][k-1].rho*mixture[i][j][k-1].u*(1-3*gas_const*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i][j][k+1].rho*mixture[i][j][k+1].u*(1-3*gas_const*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].u), dz) ;

                    double q_diff [3] = {0.0, 0.0, 0.0};
                    double q_corr [3] = {0.0, 0.0, 0.0};

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
    fill_FPC();

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
}

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
                    }
                    if(mixture[i_pb][j][k].type == TYPE_A || mixture[i_pb][j][k].type == TYPE_FS){
                        mixture[i_pb][j][k].rho = mixture[i][j][k].rho;
                        mixture[i_pb][j][k].u   = mixture[i][j][k].u;
                        mixture[i_pb][j][k].v   = mixture[i][j][k].v;
                        mixture[i_pb][j][k].w   = mixture[i][j][k].w;
                        mixture[i_pb][j][k].temp= mixture[i][j][k].temp;
                    }
                    if(mixture[i][j_nb][k].type == TYPE_A || mixture[i][j_nb][k].type == TYPE_FS){
                        mixture[i][j_nb][k].rho = mixture[i][j][k].rho;
                        mixture[i][j_nb][k].u   = mixture[i][j][k].u;
                        mixture[i][j_nb][k].v   = mixture[i][j][k].v;
                        mixture[i][j_nb][k].w   = mixture[i][j][k].w;
                        mixture[i][j_nb][k].temp= mixture[i][j][k].temp;
                    }
                    if(mixture[i][j_pb][k].type == TYPE_A || mixture[i][j_pb][k].type == TYPE_FS){
                        mixture[i][j_pb][k].rho = mixture[i][j][k].rho;
                        mixture[i][j_pb][k].u   = mixture[i][j][k].u;
                        mixture[i][j_pb][k].v   = mixture[i][j][k].v;
                        mixture[i][j_pb][k].w   = mixture[i][j][k].w;
                        mixture[i][j_pb][k].temp= mixture[i][j][k].temp;
                    }
                    if(mixture[i][j][k_nb].type == TYPE_A || mixture[i][j][k_nb].type == TYPE_FS){
                        mixture[i][j][k_nb].rho = mixture[i][j][k].rho;
                        mixture[i][j][k_nb].u   = mixture[i][j][k].u;
                        mixture[i][j][k_nb].v   = mixture[i][j][k].v;
                        mixture[i][j][k_nb].w   = mixture[i][j][k].w;
                        mixture[i][j][k_nb].temp= mixture[i][j][k].temp;
                    }
                    if(mixture[i][j][k_pb].type == TYPE_A || mixture[i][j][k_pb].type == TYPE_FS){
                        mixture[i][j][k_pb].rho = mixture[i][j][k].rho;
                        mixture[i][j][k_pb].u   = mixture[i][j][k].u;
                        mixture[i][j][k_pb].v   = mixture[i][j][k].v;
                        mixture[i][j][k_pb].w   = mixture[i][j][k].w;
                        mixture[i][j][k_pb].temp= mixture[i][j][k].temp;
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


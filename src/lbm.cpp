
#include "./headers/lbm.hpp"
#include <omp.h>

LBM::LBM(int Nx, int Ny, int Nz, double nu)
{
    this->Nx = Nx + 2;
    this->Ny = Ny + 2;
    this->Nz = Nz + 2;
    this->nu = nu;
    tau =  0.5 + 3 * nu;
    std::cout << "tau       : " << tau << std::endl;
    omega = 1.0/tau;
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
}


void LBM::Init()
{
    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (fluid1[i][j][k].type == TYPE_F || fluid1[i][j][k].type == TYPE_E)     
                {    
                    #ifdef LBM_ENTROPY
                        double zeta = GAS_CONST * TEMP;
                        for (int l = 0; l < npop; ++l)
                        {
                            double feq = fluid1[i][j][k].rho;
                                if (cx[l] == 0) feq *= (1 - (fluid1[i][j][k].u*fluid1[i][j][k].u + zeta));
                                else if (cx[l] == 1) feq *= (fluid1[i][j][k].u + (fluid1[i][j][k].u*fluid1[i][j][k].u + zeta))/2;
                                else if (cx[l] == -1) feq*= (-fluid1[i][j][k].u + (fluid1[i][j][k].u*fluid1[i][j][k].u + zeta))/2;
                            
                            #if NDIM == 2 || NDIM == 3
                                if (cy[l] == 0) feq *= (1 - (fluid1[i][j][k].v*fluid1[i][j][k].v + zeta));
                                else if (cy[l] == 1) feq *= (fluid1[i][j][k].v + (fluid1[i][j][k].v*fluid1[i][j][k].v + zeta))/2;
                                else if (cy[l] == -1) feq*= (-fluid1[i][j][k].v + (fluid1[i][j][k].v*fluid1[i][j][k].v + zeta))/2;
                            #endif

                            #if NDIM == 3
                                if (cz[l] == 0) feq *= (1 - (fluid1[i][j][k].w*fluid1[i][j][k].w + zeta));
                                else if (cz[l] == 1) feq *= (fluid1[i][j][k].w + (fluid1[i][j][k].w*fluid1[i][j][k].w + zeta))/2;
                                else if (cz[l] == -1) feq*= (-fluid1[i][j][k].w + (fluid1[i][j][k].w*fluid1[i][j][k].w + zeta))/2;
                            #endif
                            

                            fluid1[i][j][k].fpc[l]=feq;
                            fluid1[i][j][k].f[l]=feq;
                        }

                    #else
                        double uu = fluid1[i][j][k].u*fluid1[i][j][k].u + fluid1[i][j][k].v*fluid1[i][j][k].v + fluid1[i][j][k].w*fluid1[i][j][k].w;

                        for (int l = 0; l < npop; ++l)
                        {
                            double cu = cx[l]*fluid1[i][j][k].u + cy[l]*fluid1[i][j][k].v + cz[l]*fluid1[i][j][k].w;
                            double feq = wi[l]*fluid1[i][j][k].rho*(1.0+3.0*cu+4.5*cu*cu-1.5*uu);
                            fluid1[i][j][k].fpc[l]=feq;
                            fluid1[i][j][k].f[l]=feq;
                        }
                    #endif
                }
            }
        }
    }
}

void LBM::Collide_BGK()
{
    #pragma omp parallel for
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int k = 0; k < Nz; ++k)
            {
                if (fluid1[i][j][k].type==TYPE_F)
                {
                    double rho=0.0;
                    double rho_u=0.0;
                    double rho_v=0.0;
                    double rho_w=0.0;
                    //double rho_E=0.0; // NOT FINISHED YET
                                
                    for (int l = 0; l < npop; ++l)
                    {
                        rho  +=fluid1[i][j][k].f[l];
                        rho_u+=fluid1[i][j][k].f[l]*cx[l];
                        rho_v+=fluid1[i][j][k].f[l]*cy[l];
                        rho_w+=fluid1[i][j][k].f[l]*cz[l];
                    }
                    // Macroscopic quantities
                    double u=rho_u/rho;
                    double v=rho_v/rho;
                    double w=rho_w/rho;
                    
                    #ifdef LBM_ENTROPY
                        double zeta = GAS_CONST * TEMP;
                        for (int l = 0; l < npop; ++l)
                        {
                            double feq = rho;
                            double h_func = 1;
                            double B_vec[3];

                                if (cx[l] == 0) feq *= (1 - (u*u + zeta));
                                else if (cx[l] == 1) feq *= (u + (u*u + zeta))/2;
                                else if (cx[l] == -1) feq*= (-u + (u*u + zeta))/2;
                            
                            #if NDIM == 2 || NDIM == 3
                                if (cy[l] == 0) feq *= (1 - (v*v + zeta));
                                else if (cy[l] == 1) feq *= (v + (v*v + zeta))/2;
                                else if (cy[l] == -1) feq*= (-v + (v*v + zeta))/2;
                            #endif

                            #if NDIM == 3
                                if (cz[l] == 0) feq *= (1 - (w*w + zeta));
                                else if (cz[l] == 1) feq *= (w + (w*w + zeta))/2;
                                else if (cz[l] == -1) feq*= (-w + (w*w + zeta))/2;
                            #endif
                            
                            fluid1[i][j][k].fpc[l]=(1.0-omega)*fluid1[i][j][k].f[l]+omega*feq;
                        }

                    #else
                        double uu = u*u + v*v + w*w;

                        for (int l = 0; l < npop; ++l)
                        {
                            double cu=cx[l]*u+cy[l]*v+cz[l]*w;
                            double feq=wi[l]*rho*(1.0+3.0*cu+4.5*cu*cu-1.5*uu);
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
    #pragma omp parallel for
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
                        }
                        //---- Inlet/Outlet Boundary Condition
                        else if (fluid1[i_nb][j_nb][k_nb].type==TYPE_E)
                        {
                            fluid1[i][j][k].f[l] = fluid1[i_nb][j_nb][k_nb].fpc[l];
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
                        }
                    }
                }
            }
        }
    } 
}

void LBM::Quantity()
{
    #pragma omp parallel for
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
                    for (int l = 0; l < npop; ++l)
                    {
                        rho+=fluid1[i][j][k].f[l];
                        rho_u+=fluid1[i][j][k].f[l]*cx[l];
                        rho_v+=fluid1[i][j][k].f[l]*cy[l];
                        rho_w+=fluid1[i][j][k].f[l]*cz[l];
                    }
                    fluid1[i][j][k].rho=rho;
                    fluid1[i][j][k].u=rho_u/rho;
                    fluid1[i][j][k].v=rho_v/rho;
                    fluid1[i][j][k].w=rho_w/rho;
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
}
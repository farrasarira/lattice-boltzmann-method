
#include "./headers/lbm.hpp"
#include <omp.h>

LBM::LBM(int Nx, int Ny, int Nz, double nu)
{
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->nu = nu;
    tau =  0.5 + 3 * nu;
    omega = 1.0/tau;
    // allocate memory for lattice
    fluid1 = new LATTICE **[Nx];
    for (int i = 0; i < Nx; ++i)
    {
        fluid1[i] = new LATTICE *[Ny];
        for (int j = 0; j < Ny; ++j)
        {
            fluid1[i][j] = new LATTICE [Nz];
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
                #ifdef LBM_ENTROPY
                        double zeta = GAS_CONST * TEMP;
                        for (int l = 0; l < 9; ++l)
                        {
                            double feq = fluid1[i][j][k].rho;
                                if (cx[l] == 0) feq *= (1 - (fluid1[i][j][k].u*fluid1[i][j][k].u + zeta));
                                else if (cx[l] == 1) feq *= (fluid1[i][j][k].u + (fluid1[i][j][k].u*fluid1[i][j][k].u + zeta))/2;
                                else if (cx[l] == -1) feq*= (-fluid1[i][j][k].u + (fluid1[i][j][k].u*fluid1[i][j][k].u + zeta))/2;
                            
                            #if ndim == 2 || ndim == 3
                                if (cy[l] == 0) feq *= (1 - (fluid1[i][j][k].v*fluid1[i][j][k].v + zeta));
                                else if (cy[l] == 1) feq *= (fluid1[i][j][k].v + (fluid1[i][j][k].v*fluid1[i][j][k].v + zeta))/2;
                                else if (cy[l] == -1) feq*= (-fluid1[i][j][k].v + (fluid1[i][j][k].v*fluid1[i][j][k].v + zeta))/2;
                            #endif

                            #if ndim == 3
                                if (cz[l] == 0) feq *= (1 - (fluid1[i][j][k].w*fluid1[i][j][k].w + zeta));
                                else if (cz[l] == 1) feq *= (fluid1[i][j][k].w + (fluid1[i][j][k].w*fluid1[i][j][k].w + zeta))/2;
                                else if (cz[l] == -1) feq*= (-fluid1[i][j][k].w + (fluid1[i][j][k].w*fluid1[i][j][k].w + zeta))/2;
                            #endif
                            

                            fluid1[i][j][k].fpc[l]=feq;
                            fluid1[i][j][k].f[l]=feq;
                        }

                    #else
                        double uu = pow(fluid1[i][j][k].u,2)+pow(fluid1[i][j][k].v,2);

                        for (int l = 0; l < npop; ++l)
                        {
                            double cu=cx[l]*fluid1[i][j][k].u + cy[l]*fluid1[i][j][k].v;
                            double feq=w[l]*fluid1[i][j][k].rho*(1.0+3.0*cu+4.5*cu*cu-1.5*uu);
                            fluid1[i][j][k].fpc[l]=feq;
                            fluid1[i][j][k].f[l]=feq;
                        }
                    #endif
            }
        }
    }
}

void LBM::Collide_BGK()
{
    #pragma omp parallel for
    for (int i = 1; i < Nx-1; ++i)
    {
        for (int j = 1; j < Ny-1; ++j)
        {
            for (int k = 0; k < Nz; ++k)
            {
                if (fluid1[i][j][k].type==TYPE_F)
                {
                    double rho=0.0;
                    double rho_u=0.0;
                    double rho_v=0.0;
                    double rho_E=0.0; // NOT FINISHED YET
                                
                    for (int l = 0; l < npop; ++l){
                        rho+=fluid1[i][j][k].f[l];
                        rho_u+=fluid1[i][j][k].f[l]*cx[l];
                        rho_v+=fluid1[i][j][k].f[l]*cy[l];

                        rho_E+=fluid1[i][j][k].g[l];
                    }
                    // Macroscopic quantities
                    double u=rho_u/rho;
                    double v=rho_v/rho;
                    
                    #ifdef LBM_ENTROPY
                        double zeta = GAS_CONST * TEMP;
                        for (int l = 0; l < 9; ++l)
                        {
                            double feq = rho;
                            double h_func = 1;
                            double B_vec[3];

                                if (cx[l] == 0) feq *= (1 - (u*u + zeta));
                                else if (cx[l] == 1) feq *= (u + (u*u + zeta))/2;
                                else if (cx[l] == -1) feq*= (-u + (u*u + zeta))/2;
                            
                            #if ndim == 2 || ndim == 3
                                if (cy[l] == 0) feq *= (1 - (v*v + zeta));
                                else if (cy[l] == 1) feq *= (v + (v*v + zeta))/2;
                                else if (cy[l] == -1) feq*= (-v + (v*v + zeta))/2;
                            #endif

                            #if ndim == 3
                                if (cz[l] == 0) feq *= (1 - (w*w + zeta));
                                else if (cz[l] == 1) feq *= (w + (w*w + zeta))/2;
                                else if (cz[l] == -1) feq*= (-w + (w*w + zeta))/2;
                            #endif
                            
                            fluid1[i][j][k].fpc[l]=(1.0-omega)*fluid1[i][j][k].f[l]+omega*feq;
                            fluid1[i][j][k].f[l]=fluid1[i][j][k].fpc[l];
                        }

                    #else
                        double uu = pow(u,2)+pow(v,2);

                        for (int l = 0; l < npop; ++l)
                        {
                            double cu=cx[l]*u+cy[l]*v;
                            double feq=w[l]*rho*(1.0+3.0*cu+4.5*cu*cu-1.5*uu);
                            fluid1[i][j][k].fpc[l]=(1.0-omega)*fluid1[i][j][k].f[l]+omega*feq;
                            fluid1[i][j][k].f[l]=fluid1[i][j][k].fpc[l];
                        }
                    #endif

                }
            }
        }
    } 
}

void LBM::Streaming()
{
    #pragma omp parallel for
    for(int i=1; i<Nx-1; ++i)
    {
        for(int j=1; j<Ny-1; ++j)
        {
            for(int k = 0; k<Nz; ++k)
            {
                if(fluid1[i][j][k].type==TYPE_F) // 1~8
                {
                    fluid1[i][j][k].f[0]=fluid1[i][j][k].fpc[0];
                    if(fluid1[i-1][j  ][k].type==TYPE_F){ fluid1[i][j][k].f[1]=fluid1[i-1][j  ][k].fpc[ 1]; }else { fluid1[i][j][k].f[1]=fluid1[i][j][k].fpc[3]; }
                    if(fluid1[i  ][j-1][k].type==TYPE_F){ fluid1[i][j][k].f[2]=fluid1[i  ][j-1][k].fpc[ 2]; }else { fluid1[i][j][k].f[2]=fluid1[i][j][k].fpc[4]; }
                    if(fluid1[i+1][j  ][k].type==TYPE_F){ fluid1[i][j][k].f[3]=fluid1[i+1][j  ][k].fpc[ 3]; }else { fluid1[i][j][k].f[3]=fluid1[i][j][k].fpc[1]; }
                    if(fluid1[i  ][j+1][k].type==TYPE_F){ fluid1[i][j][k].f[4]=fluid1[i  ][j+1][k].fpc[ 4]; }else { fluid1[i][j][k].f[4]=fluid1[i][j][k].fpc[2]; }
                    if(fluid1[i-1][j-1][k].type==TYPE_F){ fluid1[i][j][k].f[5]=fluid1[i-1][j-1][k].fpc[ 5]; }else { fluid1[i][j][k].f[5]=fluid1[i][j][k].fpc[7]; }
                    if(fluid1[i+1][j-1][k].type==TYPE_F){ fluid1[i][j][k].f[6]=fluid1[i+1][j-1][k].fpc[ 6]; }else { fluid1[i][j][k].f[6]=fluid1[i][j][k].fpc[8]; }
                    if(fluid1[i+1][j+1][k].type==TYPE_F){ fluid1[i][j][k].f[7]=fluid1[i+1][j+1][k].fpc[ 7]; }else { fluid1[i][j][k].f[7]=fluid1[i][j][k].fpc[5]; }
                    if(fluid1[i-1][j+1][k].type==TYPE_F){ fluid1[i][j][k].f[8]=fluid1[i-1][j+1][k].fpc[ 8]; }else { fluid1[i][j][k].f[8]=fluid1[i][j][k].fpc[6]; }
                }
            }
        }
    } 
}

void LBM::BC_Noslip()
{
    #pragma omp parallel for
    for (int i = 1; i < Nx-1; ++i)
    {
        for (int k = 0; k < Nz; ++k)
        {
            // Ymin plane
            fluid1[i][1][k].f[5]=fluid1[i][1][k].fpc[8];
            fluid1[i][1][k].f[2]=fluid1[i][1][k].fpc[4];
            fluid1[i][1][k].f[6]=fluid1[i][1][k].fpc[7];

            // Ymax plane
            fluid1[i][Ny-2][k].f[8]=fluid1[i][Ny-2][k].fpc[5];
            fluid1[i][Ny-2][k].f[4]=fluid1[i][Ny-2][k].fpc[2];
            fluid1[i][Ny-2][k].f[7]=fluid1[i][Ny-2][k].fpc[6];
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
                    for (int l = 0; l < npop; ++l)
                    {
                        rho+=fluid1[i][j][k].f[l];
                        rho_u+=fluid1[i][j][k].f[l]*cx[l];
                        rho_v+=fluid1[i][j][k].f[l]*cy[l];
                    }
                    fluid1[i][j][k].rho=rho;
                    fluid1[i][j][k].u=rho_u/rho;
                    fluid1[i][j][k].v=rho_v/rho;
                }
                else
                {
                    fluid1[i][j][k].rho=1.0;
                    fluid1[i][j][k].u=0.0;
                    fluid1[i][j][k].v=0.0;
                }
            }
        }
    }
}

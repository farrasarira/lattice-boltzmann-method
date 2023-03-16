
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
    #ifdef LBM_ENTROPY

    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (fluid1[i][j][k].type == TYPE_F || fluid1[i][j][k].type == TYPE_E)     
                {   
                    for (int l = 0; l < npop; ++l)
                    {
                        double feq = fluid1[i][j][k].rho;
                        double P = 0.0;
                        
                            double eps = fluid1[i][j][k].u;
                            if (i > Nx - 3)
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].u,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].u,3)) - 4*(fluid1[i-1][j][k].rho*fluid1[i-1][j][k].u*(1-3*fluid1[i-1][j][k].gas_const*fluid1[i-1][j][k].temp)-fluid1[i-1][j][k].rho*pow(fluid1[i-1][j][k].u,3)) + (fluid1[i-2][j][k].rho*fluid1[i-2][j][k].u*(1-3*fluid1[i-2][j][k].gas_const*fluid1[i-2][j][k].temp)-fluid1[i-2][j][k].rho*pow(fluid1[i-2][j][k].u,3)))/(2*dx);
                            else
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].u,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(-3*(fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].u,3)) + 4*(fluid1[i+1][j][k].rho*fluid1[i+1][j][k].u*(1-3*fluid1[i+1][j][k].gas_const*fluid1[i+1][j][k].temp)-fluid1[i+1][j][k].rho*pow(fluid1[i+1][j][k].u,3)) - (fluid1[i+2][j][k].rho*fluid1[i+2][j][k].u*(1-3*fluid1[i+2][j][k].gas_const*fluid1[i+2][j][k].temp)-fluid1[i+2][j][k].rho*pow(fluid1[i+2][j][k].u,3)))/(2*dx);
                            
                            if (cx[l] == 0) feq *= (1 - P);
                            else if (cx[l] == 1) feq *= (eps+P)/2;
                            else if (cx[l] == -1) feq*= (-eps+P)/2;
                        
                        #if NDIM == 2 || NDIM == 3
                            eps = fluid1[i][j][k].v;
                            if (j > Ny - 3)
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].v,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].v,3)) - 4*(fluid1[i][j-1][k].rho*fluid1[i][j-1][k].v*(1-3*fluid1[i][j-1][k].gas_const*fluid1[i][j-1][k].temp)-fluid1[i][j-1][k].rho*pow(fluid1[i][j-1][k].v,3)) + (fluid1[i][j-2][k].rho*fluid1[i][j-2][k].v*(1-3*fluid1[i][j-2][k].gas_const*fluid1[i][j-2][k].temp)-fluid1[i][j-2][k].rho*pow(fluid1[i][j-2][k].v,3)))/(2*dy);
                            else
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].v,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(-3*(fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].v,3)) + 4*(fluid1[i][j+1][k].rho*fluid1[i][j+1][k].v*(1-3*fluid1[i][j+1][k].gas_const*fluid1[i][j+1][k].temp)-fluid1[i][j+1][k].rho*pow(fluid1[i][j+1][k].v,3)) - (fluid1[i][j+2][k].rho*fluid1[i][j+2][k].v*(1-3*fluid1[i][j+2][k].gas_const*fluid1[i][j+2][k].temp)-fluid1[i][j+2][k].rho*pow(fluid1[i][j+2][k].v,3)))/(2*dy);
                            
                            if (cy[l] == 0) feq *= (1 - P);
                            else if (cy[l] == 1) feq *= (eps+P)/2;
                            else if (cy[l] == -1) feq*= (-eps+P)/2;
                        #endif

                        #if NDIM == 3
                            eps = fluid1[i][j][k].w;
                            if (k > Nz - 3)
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].w,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].w,3)) - 4*(fluid1[i][j][k-1].rho*fluid1[i][j][k-1].w*(1-3*fluid1[i][j][k-1].gas_const*fluid1[i][j][k-1].temp)-fluid1[i][j][k-1].rho*pow(fluid1[i][j][k-1].w,3)) + (fluid1[i][j][k-2].rho*fluid1[i][j][k-2].w*(1-3*fluid1[i][j][k-2].gas_const*fluid1[i][j][k-2].temp)-fluid1[i][j][k-2].rho*pow(fluid1[i][j][k-2].w,3)))/(2*dz);
                            else
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].w,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega) * (-3*(fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].w,3)) + 4*(fluid1[i][j][k+1].rho*fluid1[i][j][k+1].w*(1-3*fluid1[i][j][k+1].gas_const*fluid1[i][j][k+1].temp)-fluid1[i][j][k+1].rho*pow(fluid1[i][j][k+1].w,3)) - (fluid1[i][j][k+2].rho*fluid1[i][j][k+2].w*(1-3*fluid1[i][j][k+2].gas_const*fluid1[i][j][k+2].temp)-fluid1[i][j][k+2].rho*pow(fluid1[i][j][k+2].w,3)))/(2*dz);
                            
                            if (cz[l] == 0) feq *= (1 - P);
                            else if (cz[l] == 1) feq *= (eps+P)/2;
                            else if (cz[l] == -1) feq*= (-eps+P)/2;
                        #endif
                        
                        fluid1[i][j][k].fpc[l]=feq;
                        fluid1[i][j][k].f[l]=feq;
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
    #pragma omp parallel for
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int k = 0; k < Nz; ++k)
            {
                if (fluid1[i][j][k].type==TYPE_F)
                {                   
                    #ifdef LBM_ENTROPY
                        for (int l = 0; l < npop; ++l)
                    {
                        double feq = fluid1[i][j][k].rho;
                        double P = 0.0;
                        
                            double eps = fluid1[i][j][k].u;
                            if (i > Nx - 3)
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].u,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].u,3)) - 4*(fluid1[i-1][j][k].rho*fluid1[i-1][j][k].u*(1-3*fluid1[i-1][j][k].gas_const*fluid1[i-1][j][k].temp)-fluid1[i-1][j][k].rho*pow(fluid1[i-1][j][k].u,3)) + (fluid1[i-2][j][k].rho*fluid1[i-2][j][k].u*(1-3*fluid1[i-2][j][k].gas_const*fluid1[i-2][j][k].temp)-fluid1[i-2][j][k].rho*pow(fluid1[i-2][j][k].u,3)))/(2*dx);
                            else
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].u,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(-3*(fluid1[i][j][k].rho*fluid1[i][j][k].u*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].u,3)) + 4*(fluid1[i+1][j][k].rho*fluid1[i+1][j][k].u*(1-3*fluid1[i+1][j][k].gas_const*fluid1[i+1][j][k].temp)-fluid1[i+1][j][k].rho*pow(fluid1[i+1][j][k].u,3)) - (fluid1[i+2][j][k].rho*fluid1[i+2][j][k].u*(1-3*fluid1[i+2][j][k].gas_const*fluid1[i+2][j][k].temp)-fluid1[i+2][j][k].rho*pow(fluid1[i+2][j][k].u,3)))/(2*dx);
                            
                            if (cx[l] == 0) feq *= (1 - P);
                            else if (cx[l] == 1) feq *= (eps+P)/2;
                            else if (cx[l] == -1) feq*= (-eps+P)/2;
                        
                        #if NDIM == 2 || NDIM == 3
                            eps = fluid1[i][j][k].v;
                            if (j > Ny - 3)
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].v,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].v,3)) - 4*(fluid1[i][j-1][k].rho*fluid1[i][j-1][k].v*(1-3*fluid1[i][j-1][k].gas_const*fluid1[i][j-1][k].temp)-fluid1[i][j-1][k].rho*pow(fluid1[i][j-1][k].v,3)) + (fluid1[i][j-2][k].rho*fluid1[i][j-2][k].v*(1-3*fluid1[i][j-2][k].gas_const*fluid1[i][j-2][k].temp)-fluid1[i][j-2][k].rho*pow(fluid1[i][j-2][k].v,3)))/(2*dy);
                            else
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].v,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(-3*(fluid1[i][j][k].rho*fluid1[i][j][k].v*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].v,3)) + 4*(fluid1[i][j+1][k].rho*fluid1[i][j+1][k].v*(1-3*fluid1[i][j+1][k].gas_const*fluid1[i][j+1][k].temp)-fluid1[i][j+1][k].rho*pow(fluid1[i][j+1][k].v,3)) - (fluid1[i][j+2][k].rho*fluid1[i][j+2][k].v*(1-3*fluid1[i][j+2][k].gas_const*fluid1[i][j+2][k].temp)-fluid1[i][j+2][k].rho*pow(fluid1[i][j+2][k].v,3)))/(2*dy);
                            
                            if (cy[l] == 0) feq *= (1 - P);
                            else if (cy[l] == 1) feq *= (eps+P)/2;
                            else if (cy[l] == -1) feq*= (-eps+P)/2;
                        #endif

                        #if NDIM == 3
                            eps = fluid1[i][j][k].w;
                            if (k > Nz - 3)
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].w,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega)*(3*(fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].w,3)) - 4*(fluid1[i][j][k-1].rho*fluid1[i][j][k-1].w*(1-3*fluid1[i][j][k-1].gas_const*fluid1[i][j][k-1].temp)-fluid1[i][j][k-1].rho*pow(fluid1[i][j][k-1].w,3)) + (fluid1[i][j][k-2].rho*fluid1[i][j][k-2].w*(1-3*fluid1[i][j][k-2].gas_const*fluid1[i][j][k-2].temp)-fluid1[i][j][k-2].rho*pow(fluid1[i][j][k-2].w,3)))/(2*dz);
                            else
                                P = fluid1[i][j][k].gas_const*fluid1[i][j][k].temp + pow(fluid1[i][j][k].w,2) + dt*(2-omega)/(2*fluid1[i][j][k].rho*omega) * (-3*(fluid1[i][j][k].rho*fluid1[i][j][k].w*(1-3*fluid1[i][j][k].gas_const*fluid1[i][j][k].temp)-fluid1[i][j][k].rho*pow(fluid1[i][j][k].w,3)) + 4*(fluid1[i][j][k+1].rho*fluid1[i][j][k+1].w*(1-3*fluid1[i][j][k+1].gas_const*fluid1[i][j][k+1].temp)-fluid1[i][j][k+1].rho*pow(fluid1[i][j][k+1].w,3)) - (fluid1[i][j][k+2].rho*fluid1[i][j][k+2].w*(1-3*fluid1[i][j][k+2].gas_const*fluid1[i][j][k+2].temp)-fluid1[i][j][k+2].rho*pow(fluid1[i][j][k+2].w,3)))/(2*dz);
                            
                            if (cz[l] == 0) feq *= (1 - P);
                            else if (cz[l] == 1) feq *= (eps+P)/2;
                            else if (cz[l] == -1) feq*= (-eps+P)/2;
                        #endif
                        
                        fluid1[i][j][k].fpc[l]=feq;
                        fluid1[i][j][k].f[l]=feq;
                    }

                    #else
                        double uu = pow(fluid1[i][j][k].u,2) + pow(fluid1[i][j][k].v,2) + pow(fluid1[i][j][k].w,2);

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
#ifndef LBM_H
    #define LBM_H

    #include <math.h>
    #include "setup.hpp"
    #include <iostream>
 
    class LATTICE
    {
        public:
            short type = TYPE_F;
            double f[9], fpc[9];
            double rho, u, v;
    };

    class LBM
    {
        public:
            LATTICE ** fluid1;
        
        public:
            // constructor
            LBM()
            {
                // allocate memory for lattice
                fluid1 = new LATTICE *[Nx];
                for (int i = 0; i < Nx; ++i)
                {
                    fluid1[i] = new LATTICE [Ny];
                }
            }

            // initialize equilibrium
            void Init()
            {
                for(int i = 0; i < Nx ; ++i)
                {
                    for(int j = 0; j < Ny; ++j)
                    {
                        if (fluid1[i][j].type == TYPE_F)
                        {
                            fluid1[i][j].rho = 1.0;
                            fluid1[i][j].u = U;
                            fluid1[i][j].v = 0;
                        }
                        
                        double uu = pow(fluid1[i][j].u,2)+pow(fluid1[i][j].v,2);

                        for (int l = 0; l < 9; ++l)
                        {
                            double cu=cx[l]*fluid1[i][j].u+cy[l]*fluid1[i][j].v;
                            double feq=w[l]*fluid1[i][j].rho*(1.0+3.0*cu+4.5*cu*cu-1.5*uu);
                            fluid1[i][j].f[l]=feq;
                            fluid1[i][j].fpc[l]=feq;
                        }
                    }
                }
            }

            // BGK Collision
            void Collide_BGK()
            {
                for (int i = 0; i < Nx; ++i)
                {
                    for (int j = 0; j < Ny; ++j)
                    {
                        if (fluid1[i][j].type==TYPE_F)
                        {
                            double rho=0.0;
                            double rho_u=0.0;
                            double rho_v=0.0;
                            for (int l = 0; l < 9; ++l){
                                rho+=fluid1[i][j].f[l];
                                rho_u+=fluid1[i][j].f[l]*cx[l];
                                rho_v+=fluid1[i][j].f[l]*cy[l];
                            }
                            // Macroscopic quantities
                            double u=rho_u/rho;
                            double v=rho_v/rho;
                            // Relaxation
                            // Population in post-collision state: fpc
                            double uu=u*u+v*v;
                            for (int l = 0; l < 9; ++l){
                                double cu=cx[l]*u+cy[l]*v;
                                double feq=w[l]*rho*(1.0+3.0*cu+4.5*(cu*cu)-1.5*uu);
                                fluid1[i][j].fpc[l]=(1.0-omega)*fluid1[i][j].f[l]+omega*feq;
                                fluid1[i][j].f[l]=fluid1[i][j].fpc[l];
                            }
                        }
                    }
                } 
            }

            void Streaming()
            {
               	for(int i=1; i<Nx-1; ++i)
                {
                    for(int j=1; j<Ny-1; ++j)
                    {
                        if(fluid1[i][j].type==TYPE_F) // 1~8
                        {
                            fluid1[i][j].f[0]=fluid1[i][j].fpc[0];
                            if(fluid1[i-1][j  ].type==TYPE_F){ fluid1[i][j].f[1]=fluid1[i-1][j  ].fpc[ 1]; }else { fluid1[i][j].f[1]=fluid1[i][j].fpc[3]; }
                            if(fluid1[i  ][j-1].type==TYPE_F){ fluid1[i][j].f[2]=fluid1[i  ][j-1].fpc[ 2]; }else { fluid1[i][j].f[2]=fluid1[i][j].fpc[4]; }
                            if(fluid1[i+1][j  ].type==TYPE_F){ fluid1[i][j].f[3]=fluid1[i+1][j  ].fpc[ 3]; }else { fluid1[i][j].f[3]=fluid1[i][j].fpc[1]; }
                            if(fluid1[i  ][j+1].type==TYPE_F){ fluid1[i][j].f[4]=fluid1[i  ][j+1].fpc[ 4]; }else { fluid1[i][j].f[4]=fluid1[i][j].fpc[2]; }
                            if(fluid1[i-1][j-1].type==TYPE_F){ fluid1[i][j].f[5]=fluid1[i-1][j-1].fpc[ 5]; }else { fluid1[i][j].f[5]=fluid1[i][j].fpc[7]; }
                            if(fluid1[i+1][j-1].type==TYPE_F){ fluid1[i][j].f[6]=fluid1[i+1][j-1].fpc[ 6]; }else { fluid1[i][j].f[6]=fluid1[i][j].fpc[8]; }
                            if(fluid1[i+1][j+1].type==TYPE_F){ fluid1[i][j].f[7]=fluid1[i+1][j+1].fpc[ 7]; }else { fluid1[i][j].f[7]=fluid1[i][j].fpc[5]; }
                            if(fluid1[i-1][j+1].type==TYPE_F){ fluid1[i][j].f[8]=fluid1[i-1][j+1].fpc[ 8]; }else { fluid1[i][j].f[8]=fluid1[i][j].fpc[6]; }
                        }
                    }
	            } 
            }

            void BC_Noslip()
            {
                for (int i = 0; i < Nx; ++i)
                {
                    // Ymin plane
                    fluid1[i][1].f[5]=fluid1[i][1].fpc[8];
                    fluid1[i][1].f[2]=fluid1[i][1].fpc[4];
                    fluid1[i][1].f[6]=fluid1[i][1].fpc[7];

                    // Ymax plane
                    fluid1[i][Ny].f[8]=fluid1[i][Ny].fpc[5];
                    fluid1[i][Ny].f[4]=fluid1[i][Ny].fpc[2];
                    fluid1[i][Ny].f[7]=fluid1[i][Ny].fpc[6];
                }
            }

            void Quantity(){
                for (int i = 0; i < Nx; ++i)
                {
                    for (int j = 0; j < Ny; ++j)
                    {
                        if (fluid1[i][j].type==TYPE_F)
                        {
                            // Raw moment
                            double rho=0.0;
                            double rho_u=0.0;
                            double rho_v=0.0;
                            for (int l = 0; l < 9; ++l)
                            {
                                rho+=fluid1[i][j].f[l];
                                rho_u+=fluid1[i][j].f[l]*cx[l];
                                rho_v+=fluid1[i][j].f[l]*cy[l];
                            }
                            fluid1[i][j].rho=rho;
                            fluid1[i][j].u=rho_u/rho;
                            fluid1[i][j].v=rho_v/rho;
                        }
                        else
                        {
                            fluid1[i][j].rho=1.0;
                            fluid1[i][j].u=0.0;
                            fluid1[i][j].v=0.0;
                        }
                    }
                }
            }

    };

#endif
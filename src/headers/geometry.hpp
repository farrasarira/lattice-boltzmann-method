#ifndef GEOMETRY_H
    #define GEOMETRY_H
    
    #include <math.h>
    #include "setup.hpp"
    #include "lbm.hpp"
        
    class GEOMETRY
    {
        public:
            void cylinder_generator(LBM &lb)
            {
                double rad=D/2;         // radius of cylinder
                double cx=dx*(Nx/5);    // x-center of cylinder
                double cy=dx*(Ny/2);    // y-center of cylinder

                for(int i=0; i<Nx; i++){
                    for(int j=0; j<Ny; j++){
                        double rx=dx*i-dx*0.5;
                        double ry=dx*j-dx*0.5;
                        if(sqrt(pow(cx-rx,2)+pow(cy-ry,2))<=rad){
                            lb.fluid1[i][j].type=TYPE_S;
                        }
                    }
                }
            };
    };

#endif
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
            LBM();

            // initialize
            void Init();        // initialize equilibrium           
            
            // collision operator
            void Collide_BGK(); // BGK collision
            void Collide_SMD(); // Stefan-Maxwell Diffusion collision 

            void Streaming();   // stream
            void BC_Noslip();   // no slip boundary condition
            void Quantity();    // calculate quantity

    };

#endif
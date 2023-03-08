#ifndef LBM_H
    #define LBM_H

    #include <math.h>
    #include <iostream>
    #include "defines.hpp"
 

    #ifdef D2Q9
        // D2Q9 velocity sets
        const int ndim = 2;
        const int npop = 9;
        const double cx[9] = {0.0,+1.0,-1.0, 0.0, 0.0,+1.0,-1.0,+1.0,-1.0};
        const double cy[9] = {0.0, 0.0, 0.0,+1.0,-1.0,+1.0,-1.0,-1.0,+1.0};
        const double cz[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        // D2Q9 weight factor | cs^2 = 1/3
        const double wi[9]  = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};
                                                                                                                 
    #elif defined D3Q15
        const int ndim = 3;
        const int npop = 15;
        const double cx[15] = {0.0,+1.0,-1.0, 0.0, 0.0, 0.0, 0.0,+1.0,+1.0,+1.0,+1.0,-1.0,-1.0,-1.0,-1.0};
        const double cy[15] = {0.0, 0.0, 0.0,+1.0,-1.0, 0.0, 0.0,+1.0,+1.0,-1.0,-1.0,+1.0,+1.0,-1.0,-1.0};
        const double cz[15] = {0.0, 0.0, 0.0, 0.0, 0.0,+1.0,-1.0,+1.0,-1.0,+1.0,-1.0,+1.0,-1.0,+1.0,-1.0};
        const double wi[15] = {2/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9,1/72,1/72,1/72,1/72,1/72,1/72,1/72,1/72};
    #elif defined D3Q27
        const int ndim = 3;
        const int npop = 27;
        const double cx[27] = {  0.0, +1.0, -1.0,  0.0,  0.0,  0.0,  0.0, +1.0, -1.0, +1.0, -1.0,  0.0,  0.0, +1.0, -1.0, +1.0, -1.0,  0.0,  0.0,  +1.0,  -1.0,  +1.0,  -1.0,  +1.0,  -1.0,  -1.0,  +1.0};
        const double cy[27] = {  0.0,  0.0,  0.0, +1.0, -1.0,  0.0,  0.0, +1.0, -1.0,  0.0,  0.0, +1.0, -1.0, -1.0, +1.0,  0.0,  0.0, +1.0, -1.0,  +1.0,  -1.0,  +1.0,  -1.0,  -1.0,  +1.0,  +1.0,  -1.0};
        const double cz[27] = {  0.0,  0.0,  0.0,  0.0,  0.0, +1.0, -1.0,  0.0,  0.0, +1.0, -1.0, +1.0, -1.0,  0.0,  0.0, -1.0, +1.0, -1.0, +1.0,  +1.0,  -1.0,  -1.0,  +1.0,  +1.0,  -1.0,  +1.0,  -1.0};
        const double wi[27] = {8./27,2./27,2./27,2./27,2./27,2./27,2./27,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./216,1./216,1./216,1./216,1./216,1./216,1./216,1./216};
    #endif

   

    const int opposite[27] = {0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25};

    class LATTICE
    {
        public:
            short type = TYPE_F;    // type of lattice (FLUID, SOLID, ... see setup for more)
            
            // ###### Momentum Kinetic Equation Parameter ######
            double f[npop], fpc[npop];    // distribution function, distribution function post collistion  
            double rho;           // macroscopic quantity
            double u = 0.0;       // velocity in x-direction
            double v = 0.0;       // velocity in y-direction
            double w = 0.0;       // velocity in z-direction

            /* // ####### Energy Kinetic Equation Parameter #######
            double g[npop], gpc[npop];    // energy distribution function, energy distrbution function post collision
            double temp = 0.8;      // temperature
            double totalEnergy;     // total energy (E)
            double internalEnergy;  // internal energy (U) | E = U + 1/2 * rho * u^2 */

    };

    class LBM
    {
        private:
            int Nx, Ny, Nz = 1;
            double nu, tau, omega;
        
        public:
            LATTICE *** fluid1;

        public:
            // constructor
            LBM(int Nx, int Ny, int Nz, double nu);

            // initialize
            void Init();        // initialize equilibrium           
            
            // collision operator
            void Collide_BGK(); // BGK collision

            void Streaming();   // stream
            void BC_Noslip();   // no slip boundary condition
            void Quantity();    // calculate macroscopic quantity

            int getNx(){return Nx;};
            int getNy(){return Ny;};
            int getNz(){return Nz;};
            double getNu(){return nu;};
    };

#endif
#ifndef SETUP_H
    #define SETUP_H

    #include<math.h>

    // ################ Lattice Size ##################
    #define dx 1
    #define dy 1
    #define dt 1
    #define rho0 1.0            // Density
    #define T0 1.0              // Temperature high
    #define cs 1.0/sqrt(3.0)    // Lattice sound speed

    // ################# LBM Model ####################
     #define LBM_ENTROPY      // LBM model based on entropy
    // #define LBM_CONV            // Conventional LBM

    // ################ Velocity Set ##################
    //#define D1Q3 // choose D1Q3 velocity set for 1D;
    #define D2Q9 // choose D2Q9 velocity set for 2D; 
    //#define D3Q15 // choose D3Q15 velocity set for 3D;
    //#define D3Q19 // choose D3Q19 velocity set for 3D; 
    //#define D3Q27 // choose D3Q27 velocity set for 3D;

    // ################## Dimension ###################
    #if defined D1Q3
        #define DIM 1
    #elif defined D2Q9
        #define DIM 2
    #else
        #define DIM 3
    #endif


    // ######### Single/Multicomponent flow ###########
    #define SCF // choose singlecomponent flow
    //#define MCF // choose multicomponent flow

    // ############## Boundary Condition ###############
    #define TYPE_F 0 // fluid domain
    #define TYPE_S 1 // (stationary or moving) solid boundary
    #define TYPE_E 2 // equilibrium boundary (inflow/outflow)
    #define TYPE_T 3 // temperature boundary
    #define TYPE_P 4 // periodic boundary condition


    // Parameter
    #define Th 20.0     // Temperature high [K or C]
    #define Tc 0.0      // Temperature low [K or C]
    #define Pr 0.71     // Prandtl number (rasio difusitas momentum dan konduktifitas termal)
    #define Re 100.0   // Reynolds number
    #define Ra 100000.0 // Rayleigh number
    #define beta 0.021  // Độ giãn nở nhiệt [K^-1]
    #define nu 0.002    // Kinetic viscosity
    #define GAS_CONST 1
    #define TEMP 0.4


    extern const int Nx ;
    extern const int Ny ;

    extern const int D ;  

    // D2Q9 velocity sets
    extern const double cx[9] ;
    extern const double cy[9] ;

    // D2Q9 weight factor | cs^2 
    extern const double w[9] ;

    extern const double U ;
    extern const double Ma;
    // The relaxation time for flow field
    extern const double tau ;         
    extern const double omega ;

    // Time loops
    extern const double tend ;
    extern const double tout ;
#endif


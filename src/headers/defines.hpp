#ifndef DEFINES_H
    #define DEFINES_H

    #include<math.h>


    // ################ Lattice Size ##################
    #define dx 1
    #define dy 1
    #define dt 1
    #define rho0 1.0            // Density
    #define T0 1.0              // Temperature high
    #define cs 1.0/sqrt(3.0)    // Lattice sound speed

    // ################# LBM Model ####################
    // #define LBM_ENTROPY      // LBM model based on entropy
     #define LBM_CONV            // Conventional LBM

    // ################ Velocity Set ##################
    //#define D1Q3 // choose D1Q3 velocity set for 1D;
    #define D2Q9 // choose D2Q9 velocity set for 2D; 
    //#define D3Q15 // choose D3Q15 velocity set for 3D;
    //#define D3Q19 // choose D3Q19 velocity set for 3D; 
    //#define D3Q27 // choose D3Q27 velocity set for 3D;

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
    #define GAS_CONST 1
    #define TEMP 0.4

#endif
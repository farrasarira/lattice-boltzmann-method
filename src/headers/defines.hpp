#ifndef DEFINES_H
    #define DEFINES_H

    #include<math.h>

    //_______________________________________________________________________________________________________
    // Input and choose the simulation condition

    // ################ Flow Cases ####################
    // Choose one of the following Flow Cases
    //#define CYLINDER_2D
    #define TAYLOR_GREEN_3D
    //#define TAYLOR_GREEN_2D


    // ############# Flow Parameters ##################
    // input the flow parameters
    #define RE 200.0   // Reynolds number
    #define NU 0.05     // Kinetic viscosity
    #define TEMP 0.4    // Temperature
    #define rho0 1.0            // Density
    #define T0 1.0              // Temperature high
    #define P0 0.0
    #define cs 1.0/sqrt(3.0)    // Lattice sound speed
    #define GAS_CONST 0.8

    // ########## Simulation Time & Output ############
    #define NSTEP 10000   // Maximum time step, in Lattice time unit
    #define TOUT 100        // Interval of time step to save the macroscopic quantity

    // ############# Physical Quantity ################
    #define BCXM  0.0   
    #define BCXP  0.001
    #define BCYM  0.0
    #define BCYP  0.001

    // ############# Lattice Parameters ###############
    #define dx 1
    #define dy 1
    #define dz 1
    #define dt 1
    #define NX 100
    #define NY 100
    #define NZ 100
       

    // ################# LBM Model ####################
    // Choose one of the following LBM Model
    #define LBM_ENTROPY      // LBM model based on entropy
    //#define LBM_CONV            // Conventional LBM






    // _____________________________________________________________________________________________________
    // Definition For code, don't change!

    // ################# Velocity Set ##################
    #if defined CYLINDER_2D || defined TAYLOR_GREEN_2D
        #define D2Q9
        #define NDIM 2
    #elif defined TAYLOR_GREEN_3D
        #define D3Q27
        #define NDIM 3
    #endif

    // ############## Boundary Condition ###############
    #define TYPE_F 0 // fluid domain
    #define TYPE_S 1 // (stationary or moving) solid boundary
    #define TYPE_E 2 // equilibrium boundary (inflow/outflow)
    #define TYPE_P 3 // periodic boundary
    #define TYPE_T 4 // temperature boundary


#endif
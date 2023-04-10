#ifndef DEFINES_H
    #define DEFINES_H

    #include<math.h>

    //_______________________________________________________________________________________________________
    // Input and choose the simulation condition

    // ################ Flow Cases ####################
    // Choose one of the following Flow Cases
    //#define CYLINDER_2D
    //#define TAYLOR_GREEN_2D
    #define TAYLOR_GREEN_3D
    //#define CHANNEL_FLOW_3D
    //#define CYLINDER_3D
    //#define VISCOSITY_TEST

    // ############# Flow Parameters ##################
    // input the flow parameters
    #define RE 200        // Reynolds number
    #define NU 0.0079  // Kinetic viscosity
    #define RHO0 1.0      // Density
    #define T_HIGH 1.0    // Temperature high
    #define cs 1.0/sqrt(3.0)    // Lattice sound speed
    #define TREF 1./3.
    #define PR 1.0    // Prandtl Number
    #define GAMMA 1.4 // Cp/Cv

    // ########## Simulation Time & Output ############
    #define NSTEP 1000   // Maximum time step, in Lattice time unit
    #define TOUT 10      // Interval of time step to save the macroscopic quantity

    // ############# Physical Quantity ################
    #define BCXM  0.0   
    #define BCXP  0.001
    #define BCYM  0.0
    #define BCYP  0.001

    // ############# Lattice Parameters ###############
    #define dx 1
    #define dy 1
    #define dz 1
    #define dt_sim 1
    #define NX 200
    #define NY 200
    #define NZ 200
       

    // ################# LBM Model ####################
    // Choose one of the following LBM Model
    #define LBM_EXTEND      // LBM model based on entropy
    //#define LBM_CONV            // Conventional LBM






    // _____________________________________________________________________________________________________
    // Definition For code, don't change!

    // ################# Velocity Set ##################
    #if defined CYLINDER_2D || defined TAYLOR_GREEN_2D 
        #define D2Q9
        #define NDIM 2
    #elif defined TAYLOR_GREEN_3D || defined CHANNEL_FLOW_3D || defined CYLINDER_3D || defined VISCOSITY_TEST
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
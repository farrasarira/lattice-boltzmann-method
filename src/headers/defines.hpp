#ifndef DEFINES_H
    #define DEFINES_H

    #include<math.h>

    //_______________________________________________________________________________________________________
    // Input and choose the simulation condition

    // ################ Flow Cases ####################
    // Choose one of the following Flow Cases
    //#define CYLINDER_2D
    //#define TAYLOR_GREEN_2D
    //#define TAYLOR_GREEN_3D
    //#define CHANNEL_FLOW_3D
    //#define CYLINDER_3D
    //#define VISCOSITY_TEST
    //#define SOD_SHOCK
    #define TERNARY_DIFFUSION

    // ############# Flow Parameters ##################
    // Reference Moment Value [in Lattice Unit]
    #define VEL0 0.1    // velocity    
    #define RHO0 1.0    // Density
    #define TEMP0 0.3      // Temperature

    // _____________________________________________________________________________________________________
    // Definition For code, don't change!

    // ################# Velocity Set ##################
    #if defined CYLINDER_2D || defined TAYLOR_GREEN_2D 
        #define D2Q9
        #define NDIM 2
    #elif defined TAYLOR_GREEN_3D || defined CHANNEL_FLOW_3D || defined CYLINDER_3D || defined VISCOSITY_TEST || defined SOD_SHOCK || defined TERNARY_DIFFUSION
        #define D3Q27
        #define NDIM 3
    #endif

    #if defined TERNARY_DIFFUSION
        #define MULTICOMP
    #endif

    // ############## Boundary Condition ###############
    #define TYPE_F 0 // fluid domain
    #define TYPE_S 1 // (stationary or moving) solid boundary
    #define TYPE_E 2 // equilibrium boundary (inflow/outflow)
    #define TYPE_P 3 // periodic boundary
    #define TYPE_T 4 // temperature boundary


#endif
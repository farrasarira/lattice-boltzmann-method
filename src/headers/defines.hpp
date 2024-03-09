#ifndef DEFINES_H
    #define DEFINES_H

    #include<math.h>

    // !! ALWAYS DO "MAKE CLEAN" AND THEN "MAKE" AFTER MODIFY THIS CODE !!

    //_______________________________________________________________________________________________________
    // Choose the simulation condition

    // ################ Flow Cases ####################
    // Choose one of the following Flow Cases
    // #define CYLINDER_2D
    // #define TAYLOR_GREEN_2D
    // #define TAYLOR_GREEN_3D
    // #define CHANNEL_FLOW_3D
    // #define CYLINDER_3D
    // #define VISCOSITY_TEST
    #define SOD_SHOCK
    // #define SOD_SHOCK_SIUNIT
    
    // ####### USE FD for species conservation #######
    // Uncomment for using LBM instead
    // #define FD 

    // ########## Parallel Computation ###############
    // Uncomment 2 lines of code below to utilize parallel computation using OpenMP
    #define PARALLEL
    #define NUM_THREADS 7

    // ############### OUTPUT UNIT ###################
    #define OUTPUT_SI               // uncomment for SI UNIT, comment for LATTICE UNIT

    // ############### Limiter Type ##################
    #define LIMITER_TYPE limiterVanleer
    // #define LIMITER_TYPE limiterMinmod
    // #define LIMITER_TYPE limiterMaxmod
    // #define LIMITER_TYPE limiterSuperbee
    // #define LIMITER_TYPE limiterMC


    // _____________________________________________________________________________________________________
    // Definition For code, don not change!

    // ############# Flow Parameters ##################
    // Reference Moment Value [in Lattice Unit]
    #define VEL0 0.1    // velocity    
    #define RHO0 1.0    // Density
    #define TEMP0 0.025   // Temperature

    // ################# Velocity Set ##################
    #if defined SOD_SHOCK_1D
        #define D1Q3
        #define NDIM 1    
    #elif defined TAYLOR_GREEN_2D 
        #define D2Q9
        #define NDIM 2
    #elif defined TAYLOR_GREEN_3D || defined CHANNEL_FLOW_3D || defined CYLINDER_3D || defined VISCOSITY_TEST || defined SOD_SHOCK || defined TERNARY_DIFFUSION || defined SOD_SHOCK_SIUNIT || defined SHEAR_LAYER_MULTICOMP || defined PERFECTLY_STIRRED_REACTOR_3D_AMMONIA || defined PERFECTLY_STIRRED_REACTOR_3D_HYDROGEN || defined CYLINDER_2D || defined CONDUCTION_1D || defined BPVT_1 || defined COUETTE_FLOW || defined COUETTE_FLOW_MULTICOMP || defined TAYLOR_GREEN_3D_MULTICOMP || defined COUNTERFLOW_SINGLECOMP || defined COUNTERFLOW_NONREACTIVE
        #define D3Q27
        #define NDIM 3
    #endif

    #if defined TERNARY_DIFFUSION || defined SOD_SHOCK_SIUNIT || defined SHEAR_LAYER_MULTICOMP || defined PERFECTLY_STIRRED_REACTOR_3D_AMMONIA || defined PERFECTLY_STIRRED_REACTOR_3D_HYDROGEN || defined CONDUCTION_1D || defined COUETTE_FLOW_MULTICOMP || defined TAYLOR_GREEN_3D_MULTICOMP || defined COUNTERFLOW_NONREACTIVE
        #define MULTICOMP
    #endif

    // ############## Boundary Condition ###############
    #define TYPE_F 0 // fluid domain
    #define TYPE_P 1 // periodic boundary
    #define TYPE_S 2 // dirichlet boundary condition
    #define TYPE_A 3 // Adiabatic No-Slip Wall
    #define TYPE_FS 4  // Adiabatic Free-Slip Wall
    #define TYPE_I 5 // Inflow boundary condition
    #define TYPE_O 6 // Outflow boundary condition (Neumann Boundary Condition | Zero gradient for all variables)
    

#endif
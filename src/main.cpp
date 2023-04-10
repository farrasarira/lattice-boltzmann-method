#include<iostream>

#include "./headers/lbm.hpp"
#include "./headers/setup.hpp"
#include "./headers/output.hpp"
#include<iostream>


int main()
{
    printLogo();

    // create LBM object
    LBM lb = main_setup();  
    std::cout << "-- Setup Done --" << std::endl;

    // initialize the distribution function 
    lb.Init();  
    lb.Quantity(); 
    std::cout << "-- Initialization Done --" << std::endl;

    // initialize time step & Save the macroscopic at t=0
    int step = 0;
    OutputVTK(step, lb);
    OutputKeEns(step, lb);
        
    // Simulation loop
    for (step = 1; step <= NSTEP; ++step)
    {
        lb.Collide_BGK();   // collision step
        // std::cout << "-- Collision Done --" << std::endl;
        lb.Streaming();     // streaming step & BC
        // std::cout << "-- Streaming Done --" << std::endl;
        lb.Quantity();       // Calculate macroscopic quantity
        // std::cout << "-- Calculate Quantity Done --" << std::endl;

        if (step % TOUT == 0)
        {
            //std::cout << "Step : " << step << std::endl;
            OutputKeEns(step, lb);
        }
        if (step % TOUT == 0)
        {
            OutputVTK(step, lb); // Sace the macroscopic quantity
        }
    }

    return 0;
}
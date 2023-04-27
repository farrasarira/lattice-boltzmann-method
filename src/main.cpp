#include<iostream>

#include "./headers/lbm.hpp"
#include "./headers/setup.hpp"
#include "./headers/output.hpp"
#include<iostream>


int main()
{
    printLogo();
    clock_t start;
    start = clock();

    // create LBM object
    LBM lb = main_setup();  
    std::cout << "-- Setup Done --" << std::endl;

    // initialize the distribution function 
    lb.Init(lb.mixture);  
    std::cout << "-- Initialization Done --" << std::endl;

    lb.calculate_moment(lb.mixture);

    // initialize time step & Save the macroscopic at t=0
    int step = 0;
    OutputVTK(step, lb);
    OutputKeEns(step, lb);

    // Simulation loop
    for (step = 1; step <= NSTEP; ++step)
    {
        lb.Collide(lb.mixture);   // collision step
        // std::cout << "-- Collision Done --" << std::endl;
        lb.Streaming(lb.mixture);     // streaming step & BC
        // std::cout << "-- Streaming Done --" << std::endl;

        if (step % TOUT == 0)
        {
            //std::cout << "Step : " << step << std::endl;
            OutputKeEns(step, lb);
        }
        if (step % TOUT == 0)
        {
            //OutputVTK(step, lb); // Save the macroscopic quantity
        }
    }

    std::cout << "Comp Time  : " << double(clock()-start)/double(CLOCKS_PER_SEC) << std::endl;

    return 0;
}
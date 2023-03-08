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
    std::cout << "Masuk" << std::endl;
    // initialize the distribution function 
    lb.Init();  
    std::cout << "Init Selesai" << std::endl;
    // initialize time step & Save the macroscopic at t=0
    int step = 0;
    OutputVTK(step, lb);
    OutputKeEns(lb);   
    //calcError(step, lb);    

    for (step = 1; step < NSTEP; ++step)
    {
        //std::cout << "loop start" << std::endl;
        lb.Collide_BGK();   // collision step
        //std::cout << "Collision selesai" << std::endl;
        lb.Streaming();     // streaming step & BC
        //std::cout << "Streaming selesai" << std::endl;

        if (step % TOUT == 0)
        {
            lb.Quantity();       // Calculate macroscopic quantity
            OutputVTK(step, lb); // Sace the macroscopic quantity
            std::cout << "Step : " << step << std::endl;
            OutputKeEns(lb);
            //calcError(step, lb);
        }
    }

    return 0;
}
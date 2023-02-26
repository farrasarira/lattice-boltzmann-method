#include<iostream>

#include "./headers/lbm.hpp"
#include "./headers/setup.hpp"
#include "./headers/output.hpp"
#include<iostream>


int main()
{
    printLogo();

    LBM lb = main_setup();    // create LBM object
    
    lb.Init();  // initialize the distribution function 

    int step = 0;
    int nout = 0;
    int tend = 1000000;
    int tout = 100;

    OutputVTK(step, lb);       

    for (step = 1; step < tend; ++step)
    {
        lb.Collide_BGK();   // collision step
        lb.Streaming();     // streaming step
        lb.BC_Noslip();     // no slip boundary condition
        if (step >= nout*tout)
        {
            lb.Quantity();
            OutputVTK(step, lb);
            std::cout << "Step : " << step << std::endl;
            ++nout;
        }
    }

    return 0;
}
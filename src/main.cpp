#include<iostream>

#include "./headers/lbm.hpp"
#include "./headers/geometry.hpp"
#include "./headers/output.hpp"
#include <omp.h>

int main()
{
    std::cout << R"(
             _       ____   __  __ 
            | |     |  _ \ |  \/  |
            | |     | |_) || \  / |
            | |     |  _ < | |\/| |
            | |____ | |_) || |  | |
            |______||____/ |_|  |_|
    )";
    
    std::cout << std :: endl << "      - Flow Diagnostics Laboratory ITB - " << std::endl << std::endl;
    
    LBM lb;
    cylinder_generator(lb);
    
    lb.Init();

    int step = 0;
    OutputVTK(step, lb);

    for (step = 1; step < tend; ++step)
    {
        std::cout << "Step : " << step << std::endl;
        lb.Collide_BGK();
        lb.Streaming();
        lb.BC_Noslip();
        lb.Quantity();
        OutputVTK(step, lb);
    }

    return 0;
}
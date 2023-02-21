#include<iostream>

#include"./headers/lbm.hpp"
#include"./headers/geometry.hpp"


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
    std::cout << "Memory allocated!" << std::endl;
    GEOMETRY geom;
    geom.cylinder_generator(lb);
    std::cout << "Geometry Created!" << std::endl;
    
    lb.Init();
    std::cout << "Simulation Initialized!" << std::endl;

    for (int step = 0; step < tend; ++step)
    {
        std::cout << "Step : " << step << std::endl;
        lb.Collide_BGK();
        std::cout << "Collided!" << std::endl;
        lb.Streaming();
        std::cout << "Streamed!" << std::endl;
        lb.BC_Noslip();
        std::cout << "Bounched-Back!" << std::endl;
        lb.Quantity();
        std::cout << "Quantity Updated!" << std::endl;
    }

    return 0;
}
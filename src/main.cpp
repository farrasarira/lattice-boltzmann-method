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

    // LBM
    main_setup();  

    std::cout << "Comp Time  : " << double(clock()-start)/double(CLOCKS_PER_SEC) << std::endl;

    return 0;
}
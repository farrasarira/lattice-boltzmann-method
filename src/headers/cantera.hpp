#ifndef CANTERA_H
    #define CANTERA_H

    #include "cantera/base/Solution.h"
    #include "cantera/thermo.h"
    #include "cantera/transport.h"

    #include <omp.h>

    extern int nThreads;
    extern std::vector<std::shared_ptr<Cantera::Solution>> sols;
    extern std::shared_ptr<Cantera::Solution> sol;

#endif
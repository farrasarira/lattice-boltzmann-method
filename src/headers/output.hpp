#ifndef OUTPUT_H
    #define OUTPUT_H

    #include "lbm.hpp"

    void OutputVTK(int &nout, LBM &lb);
    void printLogo();
    void OutputKeEns(LBM &lb);
    void calcError(int &t,LBM &lb);
    
#endif
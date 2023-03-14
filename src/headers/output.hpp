#ifndef OUTPUT_H
    #define OUTPUT_H

    #include "lbm.hpp"

    void OutputVTK(int &nout, LBM &lb);
    void printLogo();
    void OutputKeEns(int &, LBM &);
    void calcError(int &,LBM &);
    
#endif
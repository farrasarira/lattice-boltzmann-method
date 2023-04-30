#ifndef OUTPUT_H
    #define OUTPUT_H


    class LBM;
    void OutputVTK(int &nout, LBM *lb);
    void printLogo();
    void OutputKeEns(int &, LBM *);
    void calcError(int &,LBM *);
    
#endif
    #include "./headers/setup.hpp"
    
    extern const int Nx = 1000;
    extern const int Ny = 500;

    extern const int D = Ny/10;  

    // D2Q9 velocity sets
    extern const double cx[9] = {0.0,+1.0, 0.0,-1.0, 0.0,+1.0,-1.0,-1.0,+1.0};
    extern const double cy[9] = {0.0, 0.0,+1.0, 0.0,-1.0,+1.0,+1.0,-1.0,-1.0};

    // D2Q9 weight factor | cs^2 = 1/3
    extern const double w[9] = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};

    extern const double U = Re * nu / D;
    extern const double Ma = U * cs;
    // The relaxation time for flow field
    extern const double tau = 0.5 + 3 * nu;         
    extern const double omega = 1.0 / tau;

    // Time loops
    extern const double tend = dt * 1000000;
    extern const double tout = dt * 100;

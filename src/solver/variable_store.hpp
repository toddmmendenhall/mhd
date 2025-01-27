#pragma once

#include <vector>

namespace MHD {

struct VariableStore {
    VariableStore() = default;
    ~VariableStore() = default;

    // Constants
    double const rSpec = 8.314 / 0.0280134; // for N2
    double const gamma = 1.4;               // for N2 at 298 K and 1 atm

    // Conserved
    std::vector<double> rho;        // mass density
    std::vector<double> rho_u;      // x momentum density
    std::vector<double> rho_v;      // y momentum density
    std::vector<double> rho_w;      // z momentum density
    std::vector<double> rho_e;      // energy density
    std::vector<double> phix;       // x magnetic flux
    std::vector<double> phiy;       // y magnetic flux
    std::vector<double> phiz;       // z magnetic flux

    // Primitive
    std::vector<double> u;          // x velocity
    std::vector<double> v;          // y velocity
    std::vector<double> w;          // z velocity
    std::vector<double> eInt;       // specific internal energy
    std::vector<double> pres;       // pressure
    std::vector<double> temp;       // temperature
    std::vector<double> bx;         // x magnetic field
    std::vector<double> by;         // y magnetic field
    std::vector<double> bz;         // z magnetic field

    // Auxiliary
    std::vector<double> rhoInv;     // specific volume
    std::vector<double> u_u;        // velocity squared
    std::vector<double> b_b;        // magnetic field squared
};

} // namespace MHD
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
    std::vector<double> rhoU;       // x momentum density
    std::vector<double> rhoV;       // y momentum density
    std::vector<double> rhoW;       // z momentum density
    std::vector<double> rhoE;       // total energy density
    std::vector<double> phix;       // x magnetic flux
    std::vector<double> phiy;       // y magnetic flux
    std::vector<double> phiz;       // z magnetic flux

    // Primitive
    std::vector<double> u;          // x velocity
    std::vector<double> v;          // y velocity
    std::vector<double> w;          // z velocity
    std::vector<double> p;          // total pressure
    std::vector<double> e;          // specific internal energy
    std::vector<double> bX;         // x magnetic field
    std::vector<double> bY;         // y magnetic field
    std::vector<double> bZ;         // z magnetic field

    // Auxiliary
    std::vector<double> rhoInv;     // specific volume
    std::vector<double> uu;         // velocity squared
    std::vector<double> t;          // temperature
    std::vector<double> bb;         // magnetic field squared
    std::vector<double> eX;         // x electric field
    std::vector<double> eY;         // y electric field
    std::vector<double> eZ;         // z electric field
};

} // namespace MHD
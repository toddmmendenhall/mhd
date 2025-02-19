#pragma once

#include <vector>

namespace MHD {

struct VariableStore {
    VariableStore() = default;
    ~VariableStore() = default;

    // Constants
    double const r = 8.314 / 0.0280134; // specific gas constant for N2
    double const gamma = 1.4;           // ratio of C_p to C_v for N2 at 298 K and 1 atm

    std::size_t const numCells = 0;

    // Cell-centered
    // Conserved
    std::vector<double> rho;      // mass density
    std::vector<double> rhoU;     // x momentum density
    std::vector<double> rhoV;     // y momentum density
    std::vector<double> rhoW;     // z momentum density
    std::vector<double> rhoE;     // total energy density

    // Primitive
    std::vector<double> u;        // x velocity
    std::vector<double> v;        // y velocity
    std::vector<double> w;        // z velocity
    std::vector<double> p;        // total pressure
    std::vector<double> e;        // specific internal energy
    std::vector<double> bX;       // x magnetic field
    std::vector<double> bY;       // y magnetic field
    std::vector<double> bZ;       // z magnetic field

    // Auxiliary
    std::vector<double> rhoInv;   // specific volume
    std::vector<double> uu;       // velocity squared
    std::vector<double> t;        // temperature
    std::vector<double> bb;       // magnetic field squared

    // Face-centered
    // Primitive
    std::vector<double> faceBX;       // x magnetic field
    std::vector<double> faceBY;       // y magnetic field
    std::vector<double> faceBZ;       // z magnetic field

    // Edge-centered
    // Auxiliary
    std::vector<double> edgeEX;       // x electric field
    std::vector<double> edgeEY;       // y electric field
    std::vector<double> edgeEZ;       // z electric field
};

} // namespace MHD
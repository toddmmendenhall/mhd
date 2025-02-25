#pragma once

#include <grid.hpp>

#include <vector>

namespace MHD {

struct VariableStore {
    VariableStore(IGrid const& grid) {
        std::size_t const numCells = grid.NumCells();

        // Default conserved state
        rho.resize(numCells, 0.9);
        rhoU.resize(numCells, 0.0);
        rhoV.resize(numCells, 0.0);
        rhoW.resize(numCells, 0.0);
        rhoE.resize(numCells, 0.0);

        // Default primitve state
        u.resize(numCells, 0.0);
        v.resize(numCells, 0.0);
        w.resize(numCells, 0.0);
        p.resize(numCells, 0.0);
        e.resize(numCells, 0.0);
        bX.resize(numCells, 0.0);
        bY.resize(numCells, 0.0);
        bZ.resize(numCells, 0.0);

        // Default auxilliary variables
        rhoInv.resize(numCells, 0.0);
        uu.resize(numCells, 0.0);
        t.resize(numCells, 0.0);
        bb.resize(numCells, 0.0);

        // Default MHD variables
        faceBX.resize(numCells, 0.0);
        faceBY.resize(numCells, 0.0);
        faceBZ.resize(numCells, 0.0);
        edgeEX.resize(numCells, 0.0);
        edgeEY.resize(numCells, 0.0);
        edgeEZ.resize(numCells, 0.0);
    }
    ~VariableStore() = default;

    // Constants
    double const r = 8.314 / 0.0280134; // specific gas constant for N2 [J/(kg K)]
    double const gamma = 1.4;           // ratio of C_p to C_v for N2 at 298 K and 1 atm

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
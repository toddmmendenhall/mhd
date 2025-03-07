#pragma once

#include <constants.hpp>
#include <grid.hpp>

#include <vector>

namespace MHD {

struct VariableStore {
    VariableStore(IGrid const& grid) {
        std::size_t const size = grid.NumCells();

        // Conserved variables
        rho.resize(size, 0.0);
        rhoU.resize(size, 0.0);
        rhoV.resize(size, 0.0);
        rhoW.resize(size, 0.0);
        rhoE.resize(size, 0.0);

        // Primitve variables
        u.resize(size, 0.0);
        v.resize(size, 0.0);
        w.resize(size, 0.0);
        p.resize(size, 0.0);
        e.resize(size, 0.0);
        bX.resize(size, 0.0);
        bY.resize(size, 0.0);
        bZ.resize(size, 0.0);

        // Auxilliary variables
        rhoInv.resize(size, 0.0);
        uu.resize(size, 0.0);
        t.resize(size, 0.0);
        cs.resize(size, 0.0);
        bb.resize(size, 0.0);

        // MHD variables
        faceBX.resize(size, 0.0);
        faceBY.resize(size, 0.0);
        faceBZ.resize(size, 0.0);
        edgeEX.resize(size, 0.0);
        edgeEY.resize(size, 0.0);
        edgeEZ.resize(size, 0.0);
    }
    ~VariableStore() = default;

    // Constants
    double const r = GAS_CONSTANT / 0.0280134; // specific gas constant for N2 [J/(kg K)]
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
    std::vector<double> cs;       // sound speed squared
    std::vector<double> bb;       // magnetic field squared
    double sMax = 0.0;            // maximum wave speed

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
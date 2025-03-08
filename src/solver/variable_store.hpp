#pragma once

#include <constants.hpp>
#include <grid.hpp>

#include <vector>

namespace MHD {

struct VariableStore {
    VariableStore(IGrid const& grid) {
        std::size_t const size = grid.NumNodes();

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
        e.resize(size, 0.0);
        
        // Auxilliary variables
        p.resize(size, 0.0);
        t.resize(size, 0.0);
        cs.resize(size, 0.0);
    }
    // ~VariableStore() = default;

    // Constants
    double const r = GAS_CONSTANT / 0.0280134; // specific gas constant for N2 [J/(kg K)]
    double const gamma = 1.4;           // ratio of C_p to C_v for N2 at 298 K and 1 atm
    double sMax = 0.0;

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
    std::vector<double> e;        // specific internal energy
    
    // Auxiliary
    std::vector<double> p;        // pressure
    std::vector<double> t;        // temperature
    std::vector<double> cs;       // sound speed squared
};

} // namespace MHD
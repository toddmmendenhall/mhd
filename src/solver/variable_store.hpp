#pragma once

#include <vector>

namespace MHD {

struct VariableStore {
    VariableStore() = default;
    ~VariableStore() = default;

    // Constants
    double const r = 8.314 / 0.0280134; // specific gas constant for N2
    double const gamma = 1.4;           // ratio of C_p to C_v for N2 at 298 K and 1 atm

    // Cell-centered
    // Conserved
    std::vector<double> m_rho;      // mass density
    std::vector<double> m_rhoU;     // x momentum density
    std::vector<double> m_rhoV;     // y momentum density
    std::vector<double> m_rhoW;     // z momentum density
    std::vector<double> m_rhoE;     // total energy density

    // Primitive
    std::vector<double> m_u;        // x velocity
    std::vector<double> m_v;        // y velocity
    std::vector<double> m_w;        // z velocity
    std::vector<double> m_p;        // total pressure
    std::vector<double> m_e;        // specific internal energy
    std::vector<double> m_bX;       // x magnetic field
    std::vector<double> m_bY;       // y magnetic field
    std::vector<double> m_bZ;       // z magnetic field

    // Auxiliary
    std::vector<double> m_rhoInv;   // specific volume
    std::vector<double> m_uu;       // velocity squared
    std::vector<double> m_t;        // temperature
    std::vector<double> m_bb;       // magnetic field squared

    // Face-centered
    // Primitive
    std::vector<double> m_faceBX;       // x magnetic field
    std::vector<double> m_faceBY;       // y magnetic field
    std::vector<double> m_faceBZ;       // z magnetic field

    // Edge-centered
    // Auxiliary
    std::vector<double> m_edgeEX;       // x electric field
    std::vector<double> m_edgeEY;       // y electric field
    std::vector<double> m_edgeEZ;       // z electric field
};

} // namespace MHD
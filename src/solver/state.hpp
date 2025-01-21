#pragma once

#include <profile.hpp>

#include <vector>

namespace MHD {

struct State {
    State() = default;
    ~State() = default;
    // State(Profile const& profile, Grid const& grid) {}

    // Conserved
    std::vector<double> m_rho;
    std::vector<double> m_rho_u;
    std::vector<double> m_rho_e0;

    // Primitive
    std::vector<double> m_u;
    std::vector<double> m_energy;
    std::vector<double> m_pressure;

    // Auxiliary
    std::vector<double> m_temperature;

    std::vector<double> m_drhodx;
};

} // namespace MHD
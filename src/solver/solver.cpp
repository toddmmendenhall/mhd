#include "grid.hpp"
#include "profile.hpp"
#include "solver.hpp"
#include <state.hpp>
#include <finite_difference.hpp>

#include <memory>

namespace MHD {

Solver::Solver(Profile const& profile) {
    profile.GetSpatialDerivativeMethod();
    profile.GetTemporalIntegrationMethod();

    m_state = std::make_unique<State>();
}

Solver::~Solver() = default;

void Solver::SolveTimeStep(Grid const& grid) {
    CentralFiniteDifference1DKernel drhodxKernel(m_state->m_rho, 1.0 / grid.GetDx(), m_state->m_drhodx);

    for (std::size_t i = 1; i < grid.GetNodePositions().size() - 1; ++i) {
        drhodxKernel(i);
    }
}
    
} // namespace MHD

#include "grid.hpp"
#include "profile.hpp"
#include "solver.hpp"
#include <state.hpp>

#include <memory>

namespace MHD {

Solver::Solver(Profile const& profile) {
    profile.GetSpatialDerivativeMethod();
    profile.GetTemporalIntegrationMethod();

    m_state = std::make_unique<State>();
}

Solver::~Solver() = default;

void Solver::SolveTimeStep(Grid const& grid) {
    
}
    
} // namespace MHD

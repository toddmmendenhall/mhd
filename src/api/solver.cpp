#include "solver.hpp"
#include "solver_impl.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

Solver::Solver(std::unique_ptr<Profile> const& profile) {
    m_solverImpl = std::make_unique<SolverImpl>(profile);
}

Solver::~Solver() {}

void Solver::SolveTimeStep(std::unique_ptr<Grid> const& grid) {
    return m_solverImpl->SolveTimeStep(grid);
}


} // namespace MHD

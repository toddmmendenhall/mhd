#include <calc.hpp>
#include <grid.hpp>
#include <grid_factory.hpp>
#include <profile.hpp>
#include <solver.hpp>

#include <memory>

namespace MHD {

Calc::Calc(Profile const& profile) : m_profile(profile) {
    GridFactory gridFactory;
    m_grid = gridFactory.CreateGrid(m_profile);
    m_solver = std::make_unique<Solver>(m_profile);
}

void Calc::Run() {
    m_solver->SolveTimeStep(*m_grid);
}

} // namespace MHD

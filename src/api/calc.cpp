#include <calc.hpp>
#include <grid.hpp>
#include <profile.hpp>
#include <solver.hpp>

#include <memory>

namespace MHD {

Calc::Calc(Profile const& profile) : m_profile(profile) {
    m_grid = grid_factory(m_profile);
    m_solver = std::make_unique<Solver>(m_profile, *m_grid);
}

void Calc::Run() {}

} // namespace MHD

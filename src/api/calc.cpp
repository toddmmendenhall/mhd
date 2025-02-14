#include <calc.hpp>
#include <grid.hpp>
#include <profile.hpp>
#include <solver.hpp>

#include <memory>

namespace MHD {

Calc::Calc(Profile const& profile) : m_profile(profile) {
    m_grid = GridFactory(m_profile);
    m_solver = solverFactory(m_profile);
}

void Calc::Run() {}

} // namespace MHD

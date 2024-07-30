#include "calc.hpp"
#include "grid.hpp"
#include "profile.hpp"

#include <iostream>
#include <memory>

namespace MHD {

Calc::Calc() {
    m_profile = std::make_unique<Profile>();
}

void Calc::Run() {
    SetupCalc();
    std::cout << "Running calc..." << std::endl;
    m_solver->SolveTimeStep(m_grid);
}

void Calc::SetupCalc() {
    m_grid = std::make_unique<Grid>(m_profile);
    m_solver = std::make_unique<Solver>(m_profile);
}

} // namespace MHD

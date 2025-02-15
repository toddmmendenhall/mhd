#include <calc.hpp>
#include <execution_controller.hpp>
#include <grid.hpp>
#include <profile.hpp>
#include <solver.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

Calc::Calc(Profile const& profile) : m_profile(profile) {
    m_executionController = std::make_unique<ExecutionController>();
    m_variableStore = std::make_unique<VariableStore>();
    m_grid = gridFactory(m_profile);
    m_solver = solverFactory(m_profile, *m_executionController, *m_variableStore);
}

void Calc::Run() {
    m_solver->PerformTimeStep();
}

} // namespace MHD

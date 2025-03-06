#include <boundary_condition/boundary_condition.hpp>
#include <execution_controller.hpp>
#include <flux/flux_scheme.hpp>
#include <grid.hpp>
#include <integration/integration.hpp>
#include <kernels.hpp>
#include <profile.hpp>
#include <reconstruction/reconstruction.hpp>
#include <solver.hpp>
#include <residual.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

Solver::Solver(Profile const& profile, ExecutionController const& execCtrl, VariableStore& varStore, IGrid const& grid) :
    m_execCtrl(execCtrl), m_varStore(varStore), m_grid(grid)
{
    m_reconstruction = reconstructionFactory(profile, varStore, grid);
    m_flux = fluxFactory(profile, grid, m_reconstruction->GetContext());
    m_residual = std::make_unique<Residual>(grid, m_flux->GetContext());
    m_integrator = integratorFactory(m_residual->GetContext(), m_varStore, timeStep);
    m_boundCon = boundaryConditionFactory(profile, grid, varStore);
}

void Solver::PerformTimeStep() {
    // Use CFL condition to determine a timestep to maintain stability
    CalculateTimeStep();

    // Apply boundary conditions to ghost cell states
    m_boundCon->Compute(m_execCtrl);

    // Compute the face-centered states
    m_reconstruction->ComputeLeftRightStates(m_execCtrl);

    // Compute the face-centered fluxes
    m_flux->ComputeInterfaceFluxes(m_execCtrl);

    // Compute the cell-centered residuals
    m_residual->ComputeResidual(m_execCtrl);

    // Integrate over the timestep to update the conserved variables
    m_integrator->Compute(m_execCtrl);

    // Update the primitives
    ComputePrimitivesFromConserved();
}

void Solver::ComputePrimitivesFromConserved() {
    std::size_t const nIntCells = m_grid.NumCells();

    SpecificVolumeKernel rhoInvKern(m_varStore);
    m_execCtrl.LaunchKernel(rhoInvKern, nIntCells);

    VelocityKernel velKern(m_varStore);
    m_execCtrl.LaunchKernel(velKern, nIntCells);

    SpecificInternalEnergyKernel eKern(m_varStore);
    m_execCtrl.LaunchKernel(eKern, nIntCells);

    CaloricallyPerfectGasPressureKernel pKern(m_varStore);
    m_execCtrl.LaunchKernel(pKern, nIntCells);

    PerfectGasTemperatureKernel tKern(m_varStore);
    m_execCtrl.LaunchKernel(tKern, nIntCells);

    CaloricallyPerfectGasSoundSpeedKernel ccKern(m_varStore);
    m_execCtrl.LaunchKernel(ccKern, nIntCells);
}

void Solver::SetupConservedState() {
    std::size_t const nIntCells = m_grid.NumCells();

    SpecificVolumeKernel rhoInvKern(m_varStore);
    m_execCtrl.LaunchKernel(rhoInvKern, nIntCells);

    VelocitySquaredKernel uuKern(m_varStore);
    m_execCtrl.LaunchKernel(uuKern, nIntCells);

    CaloricallyPerfectGasSpecificInternalEnergyKernel eKern(m_varStore);
    m_execCtrl.LaunchKernel(eKern, nIntCells);

    PerfectGasTemperatureKernel tKern(m_varStore);
    m_execCtrl.LaunchKernel(tKern, nIntCells);

    CaloricallyPerfectGasSoundSpeedKernel ccKern(m_varStore);
    m_execCtrl.LaunchKernel(ccKern, nIntCells);

    MomentumDensityKernel rhoUKern(m_varStore);
    m_execCtrl.LaunchKernel(rhoUKern, nIntCells);

    TotalEnergyDensityKernel totalEnergyDensityKern(m_varStore);
    m_execCtrl.LaunchKernel(totalEnergyDensityKern, nIntCells);
}

void Solver::CalculateTimeStep() {
    CaloricallyPerfectGasSoundSpeedKernel ccKern(m_varStore);
    m_execCtrl.LaunchKernel(ccKern, m_grid.NumCells());

    MaximumWaveSpeedKernel sMaxKern(m_varStore);
    m_execCtrl.LaunchKernel(sMaxKern, m_grid.NumCells());

    timeStep = cfl * m_grid.CellSize()[0] / m_varStore.sMax;
    m_varStore.sMax = 0.0;
}

std::unique_ptr<ISolver> solverFactory(Profile const& profile, ExecutionController const& execCtrl,
                                       VariableStore& varStore, IGrid const& grid) {
    if (profile.m_compressibleOption == CompressibleOption::COMPRESSIBLE) {
        return std::make_unique<Solver>(profile, execCtrl, varStore, grid);
    }
    return nullptr;
}

} // namespace MHD

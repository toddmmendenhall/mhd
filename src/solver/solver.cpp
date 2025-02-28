#include <boundary_condition/boundary_condition.hpp>
#include <context.hpp>
#include <execution_controller.hpp>
#include <flux/flux_scheme.hpp>
#include <grid.hpp>
#include <integration.hpp>
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
    m_reconstructionContext = std::make_unique<ReconstructionContext>(m_varStore, m_grid);
    m_reconstruction = reconstructionFactory(profile, *m_reconstructionContext);
    m_fluxContext = std::make_unique<FluxContext>(m_grid, *m_reconstructionContext);
    m_fluxScheme = fluxSchemeFactory(profile, *m_fluxContext);
    m_residualContext = std::make_unique<ResidualContext>(m_grid, *m_fluxContext);
    m_residual = std::make_unique<Residual>();
    m_integrationContext = std::make_unique<IntegrationContext>(*m_residualContext, m_varStore, timeStep);
    m_integrator = integratorFactory();
    m_boundCon = boundaryConditionFactory(profile);
    m_boundaryConditionContext = std::make_unique<BoundaryConditionContext>(m_grid, m_varStore);
}

void Solver::PerformTimeStep() {
    // Use CFL condition to determine a timestep to maintain stability
    CalculateTimeStep();

    // Apply boundary conditions to ghost cell states
    ApplyBoundaryConditions();

    // Compute the face-centered states
    ReconstructVariables();

    // Compute the face-centered fluxes
    ComputeFluxes();

    // Compute the cell-centered residuals
    m_residual->Compute(m_execCtrl, *m_residualContext);

    // Integrate over the timestep to update the conserved variables
    m_integrator->Compute(m_execCtrl, *m_integrationContext);

    // Update the primitives
    ComputePrimitivesFromConserved();
}

void Solver::ComputePrimitivesFromConserved() {
    std::size_t const nIntCells = m_grid.NumNodes();

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
    std::size_t const nIntCells = m_grid.NumNodes();

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

void Solver::ComputeFluxes() {
    m_fluxScheme->ComputeInterfaceFluxes(m_execCtrl);
}

void Solver::ReconstructVariables() {
    m_reconstruction->Compute(m_execCtrl);
}

void Solver::ApplyBoundaryConditions() {
    m_boundCon->Compute(m_execCtrl, *m_boundaryConditionContext);
}

void Solver::ComputeResiduals() {
}

void Solver::CalculateTimeStep() {
    CaloricallyPerfectGasSoundSpeedKernel ccKern(m_varStore);
    m_execCtrl.LaunchKernel(ccKern, m_grid.NumNodes());

    MaximumWaveSpeedKernel sMaxKern(m_varStore);
    m_execCtrl.LaunchKernel(sMaxKern, m_grid.NumNodes());

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

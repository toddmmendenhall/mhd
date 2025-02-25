#include <boundary_condition.hpp>
#include <context.hpp>
#include <electric_field.hpp>
#include <execution_controller.hpp>
#include <flux_scheme.hpp>
#include <grid.hpp>
#include <integration.hpp>
#include <kernels.hpp>
#include <profile.hpp>
#include <reconstruction.hpp>
#include <solver.hpp>
#include <residual.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

Solver::Solver(Profile const& profile, ExecutionController const& execCtrl, VariableStore& varStore, IGrid const& grid,
    double const tStep) :
    m_execCtrl(execCtrl), m_varStore(varStore), m_grid(grid)
{
    m_bFieldCalc = std::make_unique<MagneticFieldCalculator>();
    m_eFieldCalc = std::make_unique<ElectricFieldCalculator>();
    m_reconstructionContext = std::make_unique<ReconstructionContext>(m_varStore, m_grid);
    m_reconstruction = reconstructionFactory(profile, *m_reconstructionContext);
    m_fluxContext = std::make_unique<FluxContext>(m_grid, *m_reconstructionContext);
    m_fluxScheme = fluxSchemeFactory(profile, *m_fluxContext);
    m_residualContext = std::make_unique<ResidualContext>(m_grid, *m_fluxContext);
    m_residual = std::make_unique<Residual>();
    m_integrationContext = std::make_unique<IntegrationContext>(*m_residualContext, m_varStore, tStep);
    m_integrator = integratorFactory();
    m_boundCon = boundaryConditionFactory(profile);
}

void Solver::PerformTimeStep() {
    // Apply boundary conditions to ghost cell states
    ApplyBoundaryConditions();

    // Compute the face-centered states
    ReconstructVariables();

    // Compute the face-centered fluxes
    ComputeFluxes();

    // Compute the cell-centered residuals
    m_residual->Compute(m_execCtrl, *m_residualContext);

    // Integrate over the timestep
    m_integrator->Compute(m_execCtrl, *m_integrationContext);
}

void Solver::ComputePrimitivesFromConserved() {
    SpecificVolumeKernel rhoInvKern(m_varStore);
    m_execCtrl.LaunchKernel(rhoInvKern, m_grid.NumCells());

    VelocityKernel velKern(m_varStore);
    m_execCtrl.LaunchKernel(velKern, m_grid.NumCells());

    SpecificInternalEnergyKernel eKern(m_varStore);
    m_execCtrl.LaunchKernel(eKern, m_grid.NumCells());

    CaloricallyPerfectGasPressureKernel pKern(m_varStore);
    m_execCtrl.LaunchKernel(pKern, m_grid.NumCells());

    PerfectGasTemperatureKernel tKern(m_varStore);
    m_execCtrl.LaunchKernel(tKern, m_grid.NumCells());
}

void Solver::SetupConservedState() {
    SpecificVolumeKernel rhoInvKern(m_varStore);
    m_execCtrl.LaunchKernel(rhoInvKern, m_grid.NumCells());

    VelocitySquaredKernel uuKern(m_varStore);
    m_execCtrl.LaunchKernel(uuKern, m_grid.NumInteriorCells());

    CaloricallyPerfectGasSpecificInternalEnergyKernel eKern(m_varStore);
    m_execCtrl.LaunchKernel(eKern, m_grid.NumInteriorCells());

    PerfectGasTemperatureKernel tKern(m_varStore);
    m_execCtrl.LaunchKernel(tKern, m_grid.NumInteriorCells());

    MomentumDensityKernel rhoUKern(m_varStore);
    m_execCtrl.LaunchKernel(rhoUKern, m_grid.NumInteriorCells());

    TotalEnergyDensityKernel totalEnergyDensityKern(m_varStore);
    m_execCtrl.LaunchKernel(totalEnergyDensityKern, m_varStore.rho.size());
}

void Solver::ComputeFluxes() {
    m_fluxScheme->ComputeInterfaceFluxes(m_execCtrl);
}

void Solver::ReconstructVariables() {
    m_reconstruction->Compute(m_execCtrl);
}

void Solver::ComputeElectricFields() {
    ElectricFieldContext context;
    m_eFieldCalc->Compute(m_execCtrl, context);
}

void Solver::ComputeMagneticFields() {
    MagneticFieldContext context;
    m_bFieldCalc->Compute(m_execCtrl, context);
}

void Solver::ApplyBoundaryConditions() {
    BoundaryConditionContext context(m_grid, m_varStore);
    m_boundCon->Compute(m_execCtrl, context);
}

void Solver::ComputeResiduals() {
}

std::unique_ptr<ISolver> solverFactory(Profile const& profile, ExecutionController const& execCtrl,
                                       VariableStore& varStore, IGrid const& grid, double const tStep) {
    if (profile.m_compressibleOption == CompressibleOption::COMPRESSIBLE) {
        return std::make_unique<Solver>(profile, execCtrl, varStore, grid, tStep);
    }
    return nullptr;
}

} // namespace MHD

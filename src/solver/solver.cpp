#include <boundary_condition.hpp>
#include <context.hpp>
#include <electric_field.hpp>
#include <execution_controller.hpp>
#include <flux_scheme.hpp>
#include <grid.hpp>
#include <kernels.hpp>
#include <profile.hpp>
#include <reconstruction.hpp>
#include <solver.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

Solver::Solver(Profile const& profile, ExecutionController const& execCtrl,
               VariableStore& varStore) : m_execCtrl(execCtrl), m_varStore(varStore) {
    m_bFieldCalc = std::make_unique<MagneticFieldCalculator>();
    m_eFieldCalc = std::make_unique<ElectricFieldCalculator>();
    m_fluxScheme = fluxSchemeFactory(profile);
    m_reconstruction = reconstructionFactory(profile);
    m_boundCon = boundaryConditionFactory(profile);
}

Error Solver::PerformTimeStep(IGrid const& grid) {
    ComputePrimitivesFromConserved();
    ApplyBoundaryConditions();
    ComputeFluxes(grid);
    ComputeElectricFields();
    ComputeMagneticFields();
    ReconstructVariables();
    UpdateConservedFromPrimitives();
    return Error::SUCCESS;
}

void Solver::ComputePrimitivesFromConserved()
{
    SpecificVolumeKernel rhoInvKern(m_varStore.m_rho, m_varStore.m_rhoInv);
    m_execCtrl.LaunchKernel(rhoInvKern, m_varStore.m_rho.size());

    VelocityKernel velKern(m_varStore.m_rhoInv, m_varStore.m_rhoU, m_varStore.m_rhoV, m_varStore.m_rhoW,
                           m_varStore.m_u, m_varStore.m_v, m_varStore.m_w, m_varStore.m_uu);
    m_execCtrl.LaunchKernel(velKern, m_varStore.m_rho.size());

    MagneticFieldSquaredKernel bSquaredKern(m_varStore.m_bX, m_varStore.m_bY, m_varStore.m_bZ,
                                            m_varStore.m_bb);
    m_execCtrl.LaunchKernel(bSquaredKern, m_varStore.m_rho.size());

    SpecificInternalEnergyKernel eKern(m_varStore.m_rhoInv, m_varStore.m_rhoE, m_varStore.m_uu, m_varStore.m_bb,
                                          m_varStore.m_e);
    m_execCtrl.LaunchKernel(eKern, m_varStore.m_rho.size());

    CaloricallyPerfectGasPressureKernel pKern(m_varStore.gamma, m_varStore.m_rho, m_varStore.m_e, m_varStore.m_bb,
                                                 m_varStore.m_p);
    m_execCtrl.LaunchKernel(pKern, m_varStore.m_rho.size());

    PerfectGasTemperatureKernel tempKern(m_varStore.r, m_varStore.m_rhoInv,
                                         m_varStore.m_p, m_varStore.m_t);
    m_execCtrl.LaunchKernel(tempKern, m_varStore.m_rho.size());
}

void Solver::UpdateConservedFromPrimitives() {
    MomentumDensityKernel momentumDensityKern(m_varStore.m_rho, m_varStore.m_u, m_varStore.m_v, m_varStore.m_w,
                                              m_varStore.m_rhoU, m_varStore.m_rhoV, m_varStore.m_rhoW);
    m_execCtrl.LaunchKernel(momentumDensityKern, m_varStore.m_rho.size());

    TotalEnergyDensityKernel totalEnergyDensityKern(m_varStore.m_rho, m_varStore.m_rhoInv,
                                                    m_varStore.m_uu, m_varStore.m_e,
                                                    m_varStore.m_bb, m_varStore.m_rhoE);
    m_execCtrl.LaunchKernel(totalEnergyDensityKern, m_varStore.m_rho.size());
}

void Solver::ComputeFluxes(IGrid const& grid) {
    FluxContext context(m_varStore, grid);
    m_fluxScheme->ComputeInterfaceFluxes(m_execCtrl, context);
}

void Solver::ReconstructVariables() {
    ReconstructionContext context;
    m_reconstruction->ComputeReconstructedVariables(m_execCtrl, context);
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
    BoundaryConditionContext context;
    m_boundCon->Compute(m_execCtrl, context);
}

std::unique_ptr<ISolver> solverFactory(Profile const& profile, ExecutionController const& execCtrl,
                                       VariableStore& varStore) {
    if (profile.m_compressibleOption == CompressibleOption::COMPRESSIBLE) {
        return std::make_unique<Solver>(profile, execCtrl, varStore);
    }
    return nullptr;
}

} // namespace MHD

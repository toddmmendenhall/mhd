#include <context.hpp>
#include <execution_controller.hpp>
#include <flux_scheme.hpp>
#include <kernels.hpp>
#include <profile.hpp>
#include <reconstruction.hpp>
#include <solver.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

Solver::Solver(Profile const& profile, ExecutionController const& execCtrl,
               VariableStore& varStore) : m_execCtrl(execCtrl), m_varStore(varStore) {
    m_fluxScheme = fluxSchemeFactory(profile);
    m_reconstruction = reconstructionFactory(profile);
}

Error Solver::PerformTimeStep() {
    return Error::SUCCESS;
}

void Solver::ComputePrimitivesFromConserved()
{
    SpecificVolumeKernel rhoInvKern(m_varStore.rho, m_varStore.rhoInv);
    m_execCtrl.LaunchKernel(rhoInvKern, m_varStore.rho.size());

    VelocityKernel velKern(m_varStore.rhoInv, m_varStore.rhoU, m_varStore.rhoV, m_varStore.rhoW,
                           m_varStore.u, m_varStore.v, m_varStore.w, m_varStore.uu);
    m_execCtrl.LaunchKernel(velKern, m_varStore.rho.size());

    MagneticFieldSquaredKernel bSquaredKern(m_varStore.bX, m_varStore.bY, m_varStore.bZ,
                                            m_varStore.bb);
    m_execCtrl.LaunchKernel(bSquaredKern, m_varStore.rho.size());

    SpecificInternalEnergyKernel eKern(m_varStore.rhoInv, m_varStore.rhoE, m_varStore.uu, m_varStore.bb,
                                          m_varStore.e);
    m_execCtrl.LaunchKernel(eKern, m_varStore.rho.size());

    CaloricallyPerfectGasPressureKernel pKern(m_varStore.gamma - 1.0, m_varStore.rho, m_varStore.e, m_varStore.bb,
                                                 m_varStore.p);
    m_execCtrl.LaunchKernel(pKern, m_varStore.rho.size());

    PerfectGasTemperatureKernel tempKern(1.0 / m_varStore.rSpec, m_varStore.rhoInv,
                                         m_varStore.p, m_varStore.t);
    m_execCtrl.LaunchKernel(tempKern, m_varStore.rho.size());
}

void Solver::UpdateConservedFromPrimitives() {}

void Solver::ComputeFluxes() {
    FluxContext context;
    m_fluxScheme->ComputeInterfaceFluxes(m_execCtrl, context);
}

void Solver::ReconstructVariables() {
    ReconstructionContext context;
    m_reconstruction->ComputeReconstructedVariables(m_execCtrl, context);
}

std::unique_ptr<ISolver> solverFactory(Profile const& profile, ExecutionController const& execCtrl,
                                       VariableStore& varStore) {
    if (profile.m_compressibleOption == CompressibleOption::COMPRESSIBLE) {
        return std::make_unique<Solver>(profile, execCtrl, varStore);
    }
    return nullptr;
}

} // namespace MHD

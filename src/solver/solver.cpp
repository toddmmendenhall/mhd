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

Solver::Solver(Profile const& profile) {
    varStore = std::make_unique<VariableStore>();
    execCtrl = std::make_unique<ExecutionController>();
    m_fluxScheme = flux_scheme_factory(profile, *execCtrl);
    m_reconstruction = reconstructionFactory(profile, *execCtrl);
}

Solver::~Solver() = default;

void Solver::ComputePrimitivesFromConserved()
{
    SpecificVolumeKernel rhoInvKern(varStore->rho, varStore->rhoInv);
    execCtrl->LaunchKernel(rhoInvKern, varStore->rho.size());

    VelocityKernel velKern(varStore->rhoInv, varStore->rhoU, varStore->rhoV, varStore->rhoW,
                           varStore->u, varStore->v, varStore->w, varStore->uu);
    execCtrl->LaunchKernel(velKern, varStore->rho.size());

    MagneticFieldSquaredKernel bSquaredKern(varStore->bX, varStore->bY, varStore->bZ,
                                            varStore->bb);
    execCtrl->LaunchKernel(bSquaredKern, varStore->rho.size());

    SpecificInternalEnergyKernel eKern(varStore->rhoInv, varStore->rhoE, varStore->uu, varStore->bb,
                                          varStore->e);
    execCtrl->LaunchKernel(eKern, varStore->rho.size());

    CaloricallyPerfectGasPressureKernel pKern(varStore->gamma - 1.0, varStore->rho, varStore->e, varStore->bb,
                                                 varStore->p);
    execCtrl->LaunchKernel(pKern, varStore->rho.size());

    PerfectGasTemperatureKernel tempKern(1.0 / varStore->rSpec, varStore->rhoInv,
                                         varStore->p, varStore->t);
    execCtrl->LaunchKernel(tempKern, varStore->rho.size());
}

void Solver::UpdateConservedFromPrimitives()
{
}

void Solver::ComputeFluxes() {
    FluxContext fluxContext;
    m_fluxScheme->computeInterfaceFluxes(fluxContext);
}

void Solver::ReconstructVariables() {
    ReconstructionContext context;
    m_reconstruction->compute(context);
}

} // namespace MHD

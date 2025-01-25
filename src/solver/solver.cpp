#include <execution_controller.hpp>
#include <kernels.hpp>
#include <profile.hpp>
#include <solver.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

Solver::Solver(Profile const& profile) {
    varStore = std::make_unique<VariableStore>();
    execCtrl = std::make_unique<ExecutionController>();
}

Solver::~Solver() = default;

void Solver::ComputePrimitivesFromConserved()
{
    SpecificVolumeKernel rhoInvKern(varStore->rho,
                                    varStore->rhoInv);
    execCtrl->LaunchKernel(rhoInvKern, varStore->rho.size());

    VelocityKernel velKern(varStore->rhoInv,
                           varStore->rho_u, varStore->rho_v, varStore->rho_w,
                           varStore->u, varStore->v, varStore->w, varStore->u_u);
    execCtrl->LaunchKernel(velKern, varStore->rho.size());

    SpecificInternalEnergyKernel eIntKern(varStore->rhoInv, varStore->rho_e,
                                          varStore->u_u, varStore->eInt);
    execCtrl->LaunchKernel(eIntKern, varStore->rho.size());

    CaloricallyPerfectGasPressureKernel presKern(varStore->gamma - 1.0, varStore->rho,
                                                 varStore->eInt, varStore->pres);
    execCtrl->LaunchKernel(presKern, varStore->rho.size());

    PerfectGasTemperatureKernel tempKern(1.0 / varStore->rSpec, varStore->rhoInv,
                                         varStore->pres, varStore->temp);
    execCtrl->LaunchKernel(tempKern, varStore->rho.size());
}

void Solver::UpdateConservedFromPrimitives()
{
}

} // namespace MHD

#include <context.hpp>
#include <error.hpp>
#include <execution_controller.hpp>
#include <profile.hpp>
#include <reconstruction.hpp>

#include <cstddef>

namespace MHD {

struct LinearReconstructionKernel {
    LinearReconstructionKernel(ReconstructionContext& ctx) : ctx(ctx) {}

    void operator()(std::size_t const cellIdx, std::size_t const dimIdx) {
        std::size_t const idx = cellIdx * ctx.dimension + dimIdx;
        double const coef = ctx.timeStep / ctx.gridSize[dimIdx];
        ctx.rho[idx] += coef * (ctx.rhoFluxRight[idx] - ctx.rhoFluxLeft[idx]);
        ctx.rhoU[idx] += coef * (ctx.rhoUFluxRight[idx] - ctx.rhoUFluxLeft[idx]);
        ctx.rhoV[idx] += coef * (ctx.rhoVFluxRight[idx] - ctx.rhoVFluxLeft[idx]);
        ctx.rhoW[idx] += coef * (ctx.rhoWFluxRight[idx] - ctx.rhoWFluxLeft[idx]);
        ctx.rhoE[idx] += coef * (ctx.rhoEFluxRight[idx] - ctx.rhoEFluxLeft[idx]);
        ctx.bX[idx] += coef * (ctx.bXFluxRight[idx] - ctx.bXFluxLeft[idx]);
        ctx.bY[idx] += coef * (ctx.bYFluxRight[idx] - ctx.bYFluxLeft[idx]);
        ctx.bZ[idx] += coef * (ctx.bZFluxRight[idx] - ctx.bZFluxLeft[idx]);
    }
    ReconstructionContext& ctx;
};

class LinearReconstruction : public IReconstruction {
public:
    LinearReconstruction(ExecutionController const& execCtrl) : execCtrl(execCtrl) {}
    virtual void compute(ReconstructionContext& ctx) const {
        LinearReconstructionKernel kernel(ctx);
        execCtrl.LaunchKernel(kernel, ctx.rho.size(), ctx.dimension);
    }
private:
    ExecutionController const& execCtrl;
};

std::unique_ptr<IReconstruction> reconstructionFactory(Profile const& profile, ExecutionController const& execCtrl) {
    switch (profile.GetReconstructionOption()) {
        case ReconstructionOption::LINEAR:
            return std::make_unique<LinearReconstruction>(execCtrl);
        default:
            throw Error::INVALID_RECONSTRUCTION_OPTION;
    }
}
    
} // namespace MHD

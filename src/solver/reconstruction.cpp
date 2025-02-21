#include <context.hpp>
#include <error.hpp>
#include <execution_controller.hpp>
#include <profile.hpp>
#include <reconstruction.hpp>

#include <cstddef>

namespace MHD {

struct ConstantReconstructionKernel {
    ConstantReconstructionKernel(ReconstructionContext& context)
        : m_context(context) {}

    void operator()(std::size_t const i) {
        // Get the left boundary cell
        std::size_t iLeft = i == 0 ? m_context.numInteriorCells : i - 1;
        std::size_t iRight = i;
        m_context.rhoLeft[i] = m_context.rho[iLeft];
        m_context.rhoULeft[i] = m_context.rhoU[iLeft];
        m_context.rhoVLeft[i] = m_context.rhoV[iLeft];
        m_context.rhoWLeft[i] = m_context.rhoW[iLeft];
        m_context.rhoELeft[i] = m_context.rhoE[iLeft];
        m_context.rhoRight[i] = m_context.rho[iRight];
        m_context.rhoURight[i] = m_context.rhoU[iRight];
        m_context.rhoVRight[i] = m_context.rhoV[iRight];
        m_context.rhoWRight[i] = m_context.rhoW[iRight];
        m_context.rhoERight[i] = m_context.rhoE[iRight];
    }

    ReconstructionContext& m_context;
};

class ConstantReconstruction : public IReconstruction {
public:
    ConstantReconstruction(ReconstructionContext& context)
        : m_context(context) {}
    
    void Compute(ExecutionController const& execCtrl) {
        ConstantReconstructionKernel kernel(m_context);
        execCtrl.LaunchKernel(kernel, m_context.numCells);
    }
    
    ReconstructionContext& m_context;
};

std::unique_ptr<IReconstruction> reconstructionFactory(Profile const& profile, ReconstructionContext& context) {
    if (ReconstructionOption::CONSTANT == profile.m_reconstructionOption) {
        return std::make_unique<ConstantReconstruction>(context);
    }
    throw Error::INVALID_RECONSTRUCTION_OPTION;
}
    
} // namespace MHD

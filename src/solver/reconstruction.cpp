#include <context.hpp>
#include <error.hpp>
#include <execution_controller.hpp>
#include <profile.hpp>
#include <reconstruction.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace MHD {

struct ConstantReconstructionKernel {
    ConstantReconstructionKernel(ReconstructionContext& context)
        : m_context(context) {}

    void operator()(std::size_t const i) {
        // Get the left and right cell indices for this face
        std::size_t iLeft = m_context.faceToNodeIndices[i][0];
        std::size_t iRight = m_context.faceToNodeIndices[i][0];

        m_context.rhoLeft[i] = m_context.rho[iLeft];
        m_context.uLeft[i] = m_context.u[iLeft];
        m_context.vLeft[i] = m_context.v[iLeft];
        m_context.wLeft[i] = m_context.w[iLeft];
        m_context.pLeft[i] = m_context.p[iLeft];
        m_context.eLeft[i] = m_context.e[iLeft];

        m_context.rhoRight[i] = m_context.rho[iRight];
        m_context.uRight[i] = m_context.u[iRight];
        m_context.vRight[i] = m_context.v[iRight];
        m_context.wRight[i] = m_context.w[iRight];
        m_context.pRight[i] = m_context.p[iRight];
        m_context.eRight[i] = m_context.e[iRight];
    }

    ReconstructionContext& m_context;
};

class ConstantReconstruction : public IReconstruction {
public:
    ConstantReconstruction(ReconstructionContext& context)
        : m_context(context) {}
    
    void Compute(ExecutionController const& execCtrl) {
        ConstantReconstructionKernel kernel(m_context);
        execCtrl.LaunchKernel(kernel, m_context.numFaces);
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

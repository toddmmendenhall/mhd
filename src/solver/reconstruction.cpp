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
        // Get the left cell index for this face
        std::size_t iLeft = m_context.faceToNodeIndices[i][0];

        m_context.rhoFace[i] = m_context.rho[iLeft];
        m_context.uFace[i] = m_context.u[iLeft];
        m_context.vFace[i] = m_context.v[iLeft];
        m_context.wFace[i] = m_context.w[iLeft];
        m_context.pFace[i] = m_context.p[iLeft];
        m_context.eFace[i] = m_context.e[iLeft];
    }

    ReconstructionContext& m_context;
};

struct LinearReconstructionKernel {
    LinearReconstructionKernel(ReconstructionContext& context)
        : m_context(context) {}

    void operator()(std::size_t const i) {
        // Get the left and right cell indices for this face
        std::size_t iLeft = m_context.faceToNodeIndices[i][0];
        std::size_t iRight = m_context.faceToNodeIndices[i][1];

        m_context.rhoFace[i] = 0.5 * (m_context.rho[iLeft] + m_context.rho[iRight]);
        m_context.uFace[i] = 0.5 * (m_context.u[iLeft] + m_context.u[iRight]);
        m_context.vFace[i] = 0.5 * (m_context.v[iLeft] + m_context.v[iRight]);
        m_context.wFace[i] = 0.5 * (m_context.w[iLeft] + m_context.w[iRight]);
        m_context.pFace[i] = 0.5 * (m_context.p[iLeft] + m_context.p[iRight]);
        m_context.eFace[i] = 0.5 * (m_context.e[iLeft] + m_context.e[iRight]);
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

class LinearReconstruction : public IReconstruction {
    public:
        LinearReconstruction(ReconstructionContext& context)
            : m_context(context) {}
        
        void Compute(ExecutionController const& execCtrl) {
            LinearReconstructionKernel kernel(m_context);
            execCtrl.LaunchKernel(kernel, m_context.numFaces);
        }
        
        ReconstructionContext& m_context;
    };

std::unique_ptr<IReconstruction> reconstructionFactory(Profile const& profile, ReconstructionContext& context) {
    if (ReconstructionOption::CONSTANT == profile.m_reconstructionOption) {
        return std::make_unique<ConstantReconstruction>(context);
    }
    if (ReconstructionOption::LINEAR == profile.m_reconstructionOption) {
        return std::make_unique<LinearReconstruction>(context);
    }
    throw Error::INVALID_RECONSTRUCTION_OPTION;
}
    
} // namespace MHD

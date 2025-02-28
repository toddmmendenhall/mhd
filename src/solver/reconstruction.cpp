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
        std::size_t iRight = m_context.faceToNodeIndices[i][1];

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

struct LinearReconstructionKernel {
    LinearReconstructionKernel(ReconstructionContext& context)
        : m_context(context) {}

    void operator()(std::size_t const i) {
        // Get the left and right cell indices for this face
        std::size_t iLeftMinus1 = m_context.faceToNodeIndices[i][0];
        std::size_t iLeft = m_context.faceToNodeIndices[i][1];
        std::size_t iRight = m_context.faceToNodeIndices[i][2];
        std::size_t iRightPlus1 = m_context.faceToNodeIndices[i][3];

        m_context.rhoLeft[i] = 0.5 * (m_context.rho[iLeftMinus1] + m_context.rho[iLeft]);
        m_context.uLeft[i] = 0.5 * (m_context.u[iLeftMinus1] + m_context.u[iLeft]);
        m_context.vLeft[i] = 0.5 * (m_context.v[iLeftMinus1] + m_context.v[iLeft]);
        m_context.wLeft[i] = 0.5 * (m_context.w[iLeftMinus1] + m_context.w[iLeft]);
        m_context.pLeft[i] = 0.5 * (m_context.p[iLeftMinus1] + m_context.p[iLeft]);
        m_context.eLeft[i] = 0.5 * (m_context.e[iLeftMinus1] + m_context.e[iLeft]);
        // m_context.csLeft[i] = 0.5 * (m_context.cs[iLeftMinus1] + m_context.cs[iLeft]);
        
        m_context.rhoRight[i] = 0.5 * (m_context.rho[iRight] + m_context.rho[iRightPlus1]);
        m_context.uRight[i] = 0.5 * (m_context.u[iRight] + m_context.u[iRightPlus1]);
        m_context.vRight[i] = 0.5 * (m_context.v[iRight] + m_context.v[iRightPlus1]);
        m_context.wRight[i] = 0.5 * (m_context.w[iRight] + m_context.w[iRightPlus1]);
        m_context.pRight[i] = 0.5 * (m_context.p[iRight] + m_context.p[iRightPlus1]);
        m_context.eRight[i] = 0.5 * (m_context.e[iRight] + m_context.e[iRightPlus1]);
        // m_context.csRight[i] = 0.5 * (m_context.cs[iRight] + m_context.cs[iRightPlus1]);
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

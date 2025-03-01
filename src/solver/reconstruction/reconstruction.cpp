#include <error.hpp>
#include <execution_controller.hpp>
#include <profile.hpp>
#include <reconstruction/reconstruction.hpp>
#include <variable_store.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace MHD {

struct ConstantReconstructionKernel {
    ConstantReconstructionKernel(ReconstructionContext& context) : m_context(context) {}

    void operator()(FaceIdx const i) {
        // Get the left and right cell indices for this face
        NodeIdx iLeft = m_context.faceIdxToNodeIdxs.at(i).left;
        NodeIdx iRight = m_context.faceIdxToNodeIdxs.at(i).right;

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
    LinearReconstructionKernel(ReconstructionContext& context) : m_context(context) {}

    void operator()(FaceIdx const i) {
        // Get the left and right cell indices for this face
        NodeIdx iLeft = m_context.faceIdxToNodeIdxs.at(i).left;
        NodeIdx iRight = m_context.faceIdxToNodeIdxs.at(i).right;
        NodeIdx iLeftMinusOne = m_context.faceIdxToNodeIdxs.at(i).leftMinusOne;
        NodeIdx iRightPlusOne = m_context.faceIdxToNodeIdxs.at(i).rightPlusOne;

        m_context.rhoLeft[i] = 0.5 * (m_context.rho[iLeftMinusOne] + m_context.rho[iLeft]);
        m_context.uLeft[i] = 0.5 * (m_context.u[iLeftMinusOne] + m_context.u[iLeft]);
        m_context.vLeft[i] = 0.5 * (m_context.v[iLeftMinusOne] + m_context.v[iLeft]);
        m_context.wLeft[i] = 0.5 * (m_context.w[iLeftMinusOne] + m_context.w[iLeft]);
        m_context.pLeft[i] = 0.5 * (m_context.p[iLeftMinusOne] + m_context.p[iLeft]);
        m_context.eLeft[i] = 0.5 * (m_context.e[iLeftMinusOne] + m_context.e[iLeft]);

        m_context.rhoRight[i] = 0.5 * (m_context.rho[iRight] + m_context.rho[iRightPlusOne]);
        m_context.uRight[i] = 0.5 * (m_context.u[iRight] + m_context.u[iRightPlusOne]);
        m_context.vRight[i] = 0.5 * (m_context.v[iRight] + m_context.v[iRightPlusOne]);
        m_context.wRight[i] = 0.5 * (m_context.w[iRight] + m_context.w[iRightPlusOne]);
        m_context.pRight[i] = 0.5 * (m_context.p[iRight] + m_context.p[iRightPlusOne]);
        m_context.eRight[i] = 0.5 * (m_context.e[iRight] + m_context.e[iRightPlusOne]);
    }

    ReconstructionContext& m_context;
};

ReconstructionContext::ReconstructionContext(VariableStore const& vs, IGrid const& grid) :
    rho(vs.rho), u(vs.u), v(vs.v), w(vs.w), p(vs.p), e(vs.e), cs(vs.cs),
    faceIdxToNodeIdxs(grid.GetFaceIdxToNodeIdxs()), numFaces(grid.NumFaces()) {
    rhoLeft.resize(numFaces, 0.0);
    uLeft.resize(numFaces, 0.0);
    vLeft.resize(numFaces, 0.0);
    wLeft.resize(numFaces, 0.0);
    pLeft.resize(numFaces, 0.0);
    eLeft.resize(numFaces, 0.0);
    csLeft.resize(numFaces, 0.0);
    rhoRight.resize(numFaces, 0.0);
    uRight.resize(numFaces, 0.0);
    vRight.resize(numFaces, 0.0);
    wRight.resize(numFaces, 0.0);
    pRight.resize(numFaces, 0.0);
    eRight.resize(numFaces, 0.0);
    csRight.resize(numFaces, 0.0);
}

class ConstantReconstruction : public IReconstruction {
public:
    ConstantReconstruction(VariableStore const& varStore, IGrid const& grid) {
        m_context = std::make_unique<ReconstructionContext>(varStore, grid);
    }
    
    void Compute(ExecutionController const& execCtrl) {
        ConstantReconstructionKernel kernel(*m_context);
        execCtrl.LaunchKernel(kernel, m_context->numFaces);
    }
};

class LinearReconstruction : public IReconstruction {
public:
    LinearReconstruction(VariableStore const& varStore, IGrid const& grid) {
        m_context = std::make_unique<ReconstructionContext>(varStore, grid);
    }
    
    void Compute(ExecutionController const& execCtrl) {
        LinearReconstructionKernel kernel(*m_context);
        execCtrl.LaunchKernel(kernel, m_context->numFaces);
    }
};

std::unique_ptr<IReconstruction> reconstructionFactory(Profile const& profile, VariableStore const& varStore, IGrid const& grid) {
    if (ReconstructionOption::CONSTANT == profile.m_reconstructionOption) {
        return std::make_unique<ConstantReconstruction>(varStore, grid);
    }
    if (ReconstructionOption::LINEAR == profile.m_reconstructionOption) {
        return std::make_unique<LinearReconstruction>(varStore, grid);
    }
    throw Error::INVALID_RECONSTRUCTION_OPTION;
}
    
} // namespace MHD

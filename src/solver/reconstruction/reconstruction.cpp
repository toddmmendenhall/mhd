#include <error.hpp>
#include <execution_controller.hpp>
#include <profile.hpp>
#include <reconstruction/reconstruction.hpp>
#include <variable_store.hpp>

#include <array>
#include <cstddef>
#include <cmath>
#include <vector>

namespace MHD {

struct ConstantReconstructionKernel {
    ConstantReconstructionKernel(ReconstructionContext& context) : m_context(context) {}

    void operator()(std::size_t const i) {
        std::size_t const faceIdx = m_context.faceIdxs[i];

        // Get the left and right cell indices for this face
        std::size_t const iLeft = m_context.faceIdxToNodeIdxs.at(faceIdx)[0];
        std::size_t const iRight = m_context.faceIdxToNodeIdxs.at(faceIdx)[1];

        m_context.rhoLeft[i] = m_context.rho[iLeft];
        m_context.uLeft[i] = m_context.u[iLeft];
        m_context.vLeft[i] = m_context.v[iLeft];
        m_context.wLeft[i] = m_context.w[iLeft];
        m_context.pLeft[i] = m_context.p[iLeft];
        m_context.eLeft[i] = m_context.e[iLeft];
        m_context.csLeft[i] = m_context.cs[iLeft];

        m_context.rhoRight[i] = m_context.rho[iRight];
        m_context.uRight[i] = m_context.u[iRight];
        m_context.vRight[i] = m_context.v[iRight];
        m_context.wRight[i] = m_context.w[iRight];
        m_context.pRight[i] = m_context.p[iRight];
        m_context.eRight[i] = m_context.e[iRight];
        m_context.csRight[i] = m_context.cs[iRight];
    }

    ReconstructionContext& m_context;
};

struct LinearReconstructionKernel {
    LinearReconstructionKernel(ReconstructionContext& context) : m_context(context) {}
    
    void operator()(std::size_t const i) {
        std::size_t const faceIdx = m_context.faceIdxs[i];
        
        // Get the left and right cell indices for this face
        std::size_t const iLeft = m_context.faceIdxToNodeIdxs.at(faceIdx)[0];
        std::size_t const iRight = m_context.faceIdxToNodeIdxs.at(faceIdx)[1];
        
        m_context.rhoLeft[i] = 0.5 * (m_context.rho[iLeft] + m_context.rho[iRight]);
        m_context.uLeft[i] = 0.5 * (m_context.u[iLeft] + m_context.u[iRight]);
        m_context.vLeft[i] = 0.5 * (m_context.v[iLeft] + m_context.v[iRight]);
        m_context.wLeft[i] = 0.5 * (m_context.w[iLeft] + m_context.w[iRight]);
        m_context.pLeft[i] = 0.5 * (m_context.p[iLeft] + m_context.p[iRight]);
        m_context.eLeft[i] = 0.5 * (m_context.e[iLeft] + m_context.e[iRight]);
        m_context.csLeft[i] = 0.5 * (m_context.cs[iLeft] + m_context.cs[iRight]);
        
        m_context.rhoRight[i] = 0.5 * (m_context.rho[iLeft] + m_context.rho[iRight]);
        m_context.uRight[i] = 0.5 * (m_context.u[iLeft] + m_context.u[iRight]);
        m_context.vRight[i] = 0.5 * (m_context.v[iLeft] + m_context.v[iRight]);
        m_context.wRight[i] = 0.5 * (m_context.w[iLeft] + m_context.w[iRight]);
        m_context.pRight[i] = 0.5 * (m_context.p[iLeft] + m_context.p[iRight]);
        m_context.eRight[i] = 0.5 * (m_context.e[iLeft] + m_context.e[iRight]);
        m_context.csRight[i] = 0.5 * (m_context.cs[iLeft] + m_context.cs[iRight]);
    }
    
    ReconstructionContext& m_context;
};

double vanLeer(double const r) {
    if (std::isnan(r) || std::isinf(r)) {
        return 2.0;
    }
    return (r + std::abs(r)) / (1.0 + std::abs(r));
    // return std::max(0.0, std::max(std::min(2.0 * r, 1.0), std::min(r, 2.0)));
}

struct MUSCLReconstructionKernel {
    MUSCLReconstructionKernel(ReconstructionContext& context) : m_context(context) {}

    void operator()(std::size_t const i) {
        std::size_t const faceIdx = m_context.faceIdxs[i];

        // Get the left, right, left-1, and right+1 cell indices for this face
        std::size_t const iLeft = m_context.faceIdxToNodeIdxs.at(faceIdx)[0];
        std::size_t const iRight = m_context.faceIdxToNodeIdxs.at(faceIdx)[1];
        std::size_t const iLeftMinusOne = m_context.faceIdxToNodeIdxs.at(faceIdx)[2];
        std::size_t const iRightPlusOne = m_context.faceIdxToNodeIdxs.at(faceIdx)[3];

        double rLeftRho = (m_context.rho[iLeft] - m_context.rho[iLeftMinusOne]) / (m_context.rho[iRight] - m_context.rho[iLeft]);
        double rRightRho = (m_context.rho[iRight] - m_context.rho[iLeft]) / (m_context.rho[iRightPlusOne] - m_context.rho[iRight]);
        double rLeftU = (m_context.u[iLeft] - m_context.u[iLeftMinusOne]) / (m_context.u[iRight] - m_context.u[iLeft]);
        double rRightU = (m_context.u[iRight] - m_context.u[iLeft]) / (m_context.u[iRightPlusOne] - m_context.u[iRight]);
        double rLeftV = (m_context.v[iLeft] - m_context.v[iLeftMinusOne]) / (m_context.v[iRight] - m_context.v[iLeft]);
        double rRightV = (m_context.v[iRight] - m_context.v[iLeft]) / (m_context.v[iRightPlusOne] - m_context.v[iRight]);
        double rLeftW = (m_context.w[iLeft] - m_context.w[iLeftMinusOne]) / (m_context.w[iRight] - m_context.w[iLeft]);
        double rRightW = (m_context.w[iRight] - m_context.w[iLeft]) / (m_context.w[iRightPlusOne] - m_context.w[iRight]);
        double rLeftP = (m_context.p[iLeft] - m_context.p[iLeftMinusOne]) / (m_context.p[iRight] - m_context.p[iLeft]);
        double rRightP = (m_context.p[iRight] - m_context.p[iLeft]) / (m_context.p[iRightPlusOne] - m_context.p[iRight]);
        double rLeftE = (m_context.e[iLeft] - m_context.e[iLeftMinusOne]) / (m_context.e[iRight] - m_context.e[iLeft]);
        double rRightE = (m_context.e[iRight] - m_context.e[iLeft]) / (m_context.e[iRightPlusOne] - m_context.e[iRight]);
        double rLeftCs = (m_context.cs[iLeft] - m_context.cs[iLeftMinusOne]) / (m_context.cs[iRight] - m_context.cs[iLeft]);
        double rRightCs = (m_context.cs[iRight] - m_context.cs[iLeft]) / (m_context.cs[iRightPlusOne] - m_context.cs[iRight]);

        double phiLeftRho = vanLeer(rLeftRho);
        double phiRightRho = vanLeer(rRightRho);
        double phiLeftU = vanLeer(rLeftU);
        double phiRightU = vanLeer(rRightU);
        double phiLeftV = vanLeer(rLeftV);
        double phiRightV = vanLeer(rRightV);
        double phiLeftW = vanLeer(rLeftW);
        double phiRightW = vanLeer(rRightW);
        double phiLeftP = vanLeer(rLeftP);
        double phiRightP = vanLeer(rRightP);
        double phiLeftE = vanLeer(rLeftE);
        double phiRightE = vanLeer(rRightE);
        double phiLeftCs = vanLeer(rLeftCs);
        double phiRightCs = vanLeer(rRightCs);

        m_context.rhoLeft[i] = m_context.rho[iLeft] + 0.5 * phiLeftRho * (m_context.rho[iRight] - m_context.rho[iLeft]);
        m_context.uLeft[i] = m_context.u[iLeft] + 0.5 * phiLeftU * (m_context.u[iRight] - m_context.u[iLeft]);
        m_context.vLeft[i] = m_context.v[iLeft] + 0.5 * phiLeftV * (m_context.v[iRight] - m_context.v[iLeft]);
        m_context.wLeft[i] = m_context.w[iLeft] + 0.5 * phiLeftW * (m_context.w[iRight] - m_context.w[iLeft]);
        m_context.pLeft[i] = m_context.p[iLeft] + 0.5 * phiLeftP * (m_context.p[iRight] - m_context.p[iLeft]);
        m_context.eLeft[i] = m_context.e[iLeft] + 0.5 * phiLeftE * (m_context.e[iRight] - m_context.e[iLeft]);
        m_context.csLeft[i] = m_context.cs[iLeft] + 0.5 * phiLeftCs * (m_context.cs[iRight] - m_context.cs[iLeft]);

        m_context.rhoRight[i] = m_context.rho[iRight] - 0.5 * phiRightRho * (m_context.rho[iRightPlusOne] - m_context.rho[iRight]);
        m_context.uRight[i] = m_context.u[iRight] - 0.5 * phiRightU * (m_context.u[iRightPlusOne] - m_context.u[iRight]);
        m_context.vRight[i] = m_context.v[iRight] - 0.5 * phiRightV * (m_context.v[iRightPlusOne] - m_context.v[iRight]);
        m_context.wRight[i] = m_context.w[iRight] - 0.5 * phiRightW * (m_context.w[iRightPlusOne] - m_context.w[iRight]);
        m_context.pRight[i] = m_context.p[iRight] - 0.5 * phiRightP * (m_context.p[iRightPlusOne] - m_context.p[iRight]);
        m_context.eRight[i] = m_context.e[iRight] - 0.5 * phiRightE * (m_context.e[iRightPlusOne] - m_context.e[iRight]);
        m_context.csRight[i] = m_context.cs[iRight] - 0.5 * phiRightCs * (m_context.cs[iRightPlusOne] - m_context.cs[iRight]);
    }

    ReconstructionContext& m_context;
};

ReconstructionContext::ReconstructionContext(VariableStore const& vs, IGrid const& grid) :
    rho(vs.rho), u(vs.u), v(vs.v), w(vs.w), p(vs.p), e(vs.e), cs(vs.cs),
    faceIdxToNodeIdxs(grid.FaceIdxToCellIdxs()), numFaces(grid.NumFaces()), faceIdxs(grid.FaceIdxs()) {
        std::size_t const size = numFaces;
        rhoLeft.resize(size, 0.0);
        uLeft.resize(size, 0.0);
        vLeft.resize(size, 0.0);
        wLeft.resize(size, 0.0);
        pLeft.resize(size, 0.0);
        eLeft.resize(size, 0.0);
        csLeft.resize(size, 0.0);
        rhoRight.resize(size, 0.0);
        uRight.resize(size, 0.0);
        vRight.resize(size, 0.0);
        wRight.resize(size, 0.0);
        pRight.resize(size, 0.0);
        eRight.resize(size, 0.0);
        csRight.resize(size, 0.0);
}

class ConstantReconstruction : public IReconstruction {
public:
    ConstantReconstruction(VariableStore const& varStore, IGrid const& grid) {
        m_context = std::make_unique<ReconstructionContext>(varStore, grid);
    }
    
    void ComputeLeftRightStates(ExecutionController const& execCtrl) {
        ConstantReconstructionKernel kernel(*m_context);
        execCtrl.LaunchKernel(kernel, m_context->numFaces);
    }
};

class LinearReconstruction : public IReconstruction {
public:
    LinearReconstruction(VariableStore const& varStore, IGrid const& grid) {
        m_context = std::make_unique<ReconstructionContext>(varStore, grid);
    }
    
    void ComputeLeftRightStates(ExecutionController const& execCtrl) {
        LinearReconstructionKernel kernel(*m_context);
        execCtrl.LaunchKernel(kernel, m_context->numFaces);
    }
};

class MUSCLReconstruction : public IReconstruction {
    public:
        MUSCLReconstruction(VariableStore const& varStore, IGrid const& grid) {
            m_context = std::make_unique<ReconstructionContext>(varStore, grid);
        }
        
        void ComputeLeftRightStates(ExecutionController const& execCtrl) {
            MUSCLReconstructionKernel kernel(*m_context);
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
    if (ReconstructionOption::MUSCL == profile.m_reconstructionOption) {
        return std::make_unique<MUSCLReconstruction>(varStore, grid);
    }
    throw Error::INVALID_RECONSTRUCTION_OPTION;
}
    
} // namespace MHD

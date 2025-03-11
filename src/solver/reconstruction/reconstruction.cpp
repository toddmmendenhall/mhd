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
        std::size_t const iFace = m_context.faceIdxs.at(i);

        std::size_t const iLeft = m_context.faceIdxToNodeIdxs.at(i)[0];
        std::size_t const iRight = m_context.faceIdxToNodeIdxs.at(i)[1];
        std::size_t const iLeftMinusOne = m_context.faceIdxToNodeIdxs.at(i)[2];
        std::size_t const iRightPlusOne = m_context.faceIdxToNodeIdxs.at(i)[3];

        double rLeftRho = (m_context.rho[iRight] - m_context.rho[iLeft]) / (m_context.rho[iLeft] - m_context.rho[iLeftMinusOne]);
        double rLeftU = (m_context.u[iRight] - m_context.u[iLeft]) / (m_context.u[iLeft] - m_context.u[iLeftMinusOne]);
        double rLeftV = (m_context.v[iRight] - m_context.v[iLeft]) / (m_context.v[iLeft] - m_context.v[iLeftMinusOne]);
        double rLeftW = (m_context.w[iRight] - m_context.w[iLeft]) / (m_context.w[iLeft] - m_context.w[iLeftMinusOne]);
        double rLeftP = (m_context.p[iRight] - m_context.p[iLeft]) / (m_context.p[iLeft] - m_context.p[iLeftMinusOne]);
        double rLeftE = (m_context.e[iRight] - m_context.e[iLeft]) / (m_context.e[iLeft] - m_context.e[iLeftMinusOne]);
        double rLeftCs = (m_context.cs[iRight] - m_context.cs[iLeft]) / (m_context.cs[iLeft] - m_context.cs[iLeftMinusOne]);
        
        double rRightRho = (m_context.rho[iRightPlusOne] - m_context.rho[iRight]) / (m_context.rho[iRight] - m_context.rho[iLeft]);
        double rRightU = (m_context.u[iRightPlusOne] - m_context.u[iRight]) / (m_context.u[iRight] - m_context.u[iLeft]);
        double rRightV = (m_context.v[iRightPlusOne] - m_context.v[iRight]) / (m_context.v[iRight] - m_context.v[iLeft]);
        double rRightW = (m_context.w[iRightPlusOne] - m_context.w[iRight]) / (m_context.w[iRight] - m_context.w[iLeft]);
        double rRightP = (m_context.p[iRightPlusOne] - m_context.p[iRight]) / (m_context.p[iRight] - m_context.p[iLeft]);
        double rRightE = (m_context.e[iRightPlusOne] - m_context.e[iRight]) / (m_context.e[iRight] - m_context.e[iLeft]);
        double rRightCs = (m_context.cs[iRightPlusOne] - m_context.cs[iRight]) / (m_context.cs[iRight] - m_context.cs[iLeft]);

        double phiLeftRho = vanLeer(rLeftRho);
        double phiLeftU = vanLeer(rLeftU);
        double phiLeftV = vanLeer(rLeftV);
        double phiLeftW = vanLeer(rLeftW);
        double phiLeftP = vanLeer(rLeftP);
        double phiLeftE = vanLeer(rLeftE);
        double phiLeftCs = vanLeer(rLeftCs);

        double phiLeftRhoInv = vanLeer(1.0 / rLeftRho);
        double phiLeftUInv = vanLeer(1.0 / rLeftU);
        double phiLeftVInv = vanLeer(1.0 / rLeftV);
        double phiLeftWInv = vanLeer(1.0 / rLeftW);
        double phiLeftPInv = vanLeer(1.0 / rLeftP);
        double phiLeftEInv = vanLeer(1.0 / rLeftE);
        double phiLeftCsInv = vanLeer(1.0 / rLeftCs);
        
        double phiRightRho = vanLeer(rRightRho);
        double phiRightU = vanLeer(rRightU);
        double phiRightV = vanLeer(rRightV);
        double phiRightW = vanLeer(rRightW);
        double phiRightP = vanLeer(rRightP);
        double phiRightE = vanLeer(rRightE);
        double phiRightCs = vanLeer(rRightCs);

        double phiRightRhoInv = vanLeer(1.0 / rRightRho);
        double phiRightUInv = vanLeer(1.0 / rRightU);
        double phiRightVInv = vanLeer(1.0 / rRightV);
        double phiRightWInv = vanLeer(1.0 / rRightW);
        double phiRightPInv = vanLeer(1.0 / rRightP);
        double phiRightEInv = vanLeer(1.0 / rRightE);
        double phiRightCsInv = vanLeer(1.0 / rRightCs);

        m_context.rhoLeft[i] = m_context.rho[iLeft] + 0.25 * m_phi * ((1.0 - m_kappa) * phiLeftRho * (m_context.rho[iLeft] - m_context.rho[iLeftMinusOne]) +
                                                                      (1.0 + m_kappa) * phiLeftRhoInv * (m_context.rho[iRight] - m_context.rho[iLeft]));
        m_context.uLeft[i] = m_context.u[iLeft] + 0.25 * m_phi * ((1.0 - m_kappa) * phiLeftU * (m_context.u[iLeft] - m_context.u[iLeftMinusOne]) +
                                                                  (1.0 + m_kappa) * phiLeftUInv * (m_context.u[iRight] - m_context.u[iLeft]));
        m_context.vLeft[i] = m_context.v[iLeft] + 0.25 * m_phi * ((1.0 - m_kappa) * phiLeftV * (m_context.v[iLeft] - m_context.v[iLeftMinusOne]) +
                                                                  (1.0 + m_kappa) * phiLeftVInv * (m_context.v[iRight] - m_context.v[iLeft]));
        m_context.wLeft[i] = m_context.w[iLeft] + 0.25 * m_phi * ((1.0 - m_kappa) * phiLeftW * (m_context.w[iLeft] - m_context.w[iLeftMinusOne]) +
                                                                  (1.0 + m_kappa) * phiLeftWInv * (m_context.v[iRight] - m_context.v[iLeft]));
        m_context.pLeft[i] = m_context.p[iLeft] + 0.25 * m_phi * ((1.0 - m_kappa) * phiLeftP * (m_context.p[iLeft] - m_context.p[iLeftMinusOne]) +
                                                                  (1.0 + m_kappa) * phiLeftPInv * (m_context.p[iRight] - m_context.p[iLeft]));
        m_context.eLeft[i] = m_context.e[iLeft] + 0.25 * m_phi * ((1.0 - m_kappa) * phiLeftE * (m_context.e[iLeft] - m_context.e[iLeftMinusOne]) +
                                                                  (1.0 + m_kappa) * phiLeftEInv * (m_context.e[iRight] - m_context.e[iLeft]));
        m_context.csLeft[i] = m_context.cs[iLeft] + 0.25 * m_phi * ((1.0 - m_kappa) * phiLeftCs * (m_context.cs[iLeft] - m_context.cs[iLeftMinusOne]) +
                                                                    (1.0 + m_kappa) * phiLeftCs * (m_context.cs[iRight] - m_context.cs[iLeft]));


        m_context.rhoRight[i] = m_context.rho[iRight] - 0.25 * m_phi * ((1.0 + m_kappa) * phiRightRho * (m_context.rho[iRight] - m_context.rho[iLeft]) +
                                                                       (1.0 - m_kappa) * phiRightRhoInv * (m_context.rho[iRightPlusOne] - m_context.rho[iRight]));
        m_context.uRight[i] = m_context.u[iRight] - 0.25 * m_phi * ((1.0 + m_kappa) * phiRightU * (m_context.u[iRight] - m_context.u[iLeft]) +
                                                                   (1.0 - m_kappa) * phiRightUInv * (m_context.u[iRightPlusOne] - m_context.u[iRight]));
        m_context.vRight[i] = m_context.v[iRight] - 0.25 * m_phi * ((1.0 + m_kappa) * phiRightV * (m_context.v[iRight] - m_context.v[iLeft]) +
                                                                   (1.0 - m_kappa) * phiRightVInv * (m_context.v[iRightPlusOne] - m_context.v[iRight]));
        m_context.wRight[i] = m_context.w[iRight] - 0.25 * m_phi * ((1.0 + m_kappa) * phiRightW * (m_context.w[iRight] - m_context.w[iLeft]) +
                                                                   (1.0 - m_kappa) * phiRightWInv * (m_context.v[iRightPlusOne] - m_context.v[iRight]));
        m_context.pRight[i] = m_context.p[iRight] - 0.25 * m_phi * ((1.0 + m_kappa) * phiRightP * (m_context.p[iRight] - m_context.p[iLeft]) +
                                                                   (1.0 - m_kappa) * phiRightPInv * (m_context.p[iRightPlusOne] - m_context.p[iRight]));
        m_context.eRight[i] = m_context.e[iRight] - 0.25 * m_phi * ((1.0 + m_kappa) * phiRightE * (m_context.e[iRight] - m_context.e[iLeft]) +
                                                                   (1.0 - m_kappa) * phiRightEInv * (m_context.e[iRightPlusOne] - m_context.e[iRight]));
        m_context.csRight[i] = m_context.cs[iRight] - 0.25 * m_phi * ((1.0 + m_kappa) * phiRightCs * (m_context.cs[iRight] - m_context.cs[iLeft]) +
                                                                     (1.0 - m_kappa) * phiRightCsInv * (m_context.cs[iRightPlusOne] - m_context.cs[iRight]));
    }

    ReconstructionContext& m_context;
    double const m_phi = 1.0;
    double const m_kappa = -1.0;
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

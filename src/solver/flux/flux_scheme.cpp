#include <constants.hpp>
#include <error.hpp>
#include <execution_controller.hpp>
#include <flux/flux_scheme.hpp>
#include <grid.hpp>
#include <profile.hpp>
#include <profile_options.hpp>
#include <reconstruction/reconstruction.hpp>

#include <cstdlib>

namespace MHD {

struct KTFluxKernel {
    KTFluxKernel(FluxContext& context) : m_context(context) {}

    inline void operator()(std::size_t const i) {
        std::size_t const faceIdx = m_context.faceIdxs[i];
    
        // Face-centered magnitude of the velocity on the left
        auto uuLeft =  m_context.uLeft[faceIdx] * m_context.uLeft[faceIdx] +
                       m_context.vLeft[faceIdx] * m_context.vLeft[faceIdx] +
                       m_context.wLeft[faceIdx] * m_context.wLeft[faceIdx];

        // Face-centered magnitude of the velocity on the right
        auto uuRight =  m_context.uRight[faceIdx] * m_context.uRight[faceIdx] +
                        m_context.vRight[faceIdx] * m_context.vRight[faceIdx] +
                        m_context.wRight[faceIdx] * m_context.wRight[faceIdx];

        // Face-centered magnitude of the magnetic field on the left
        auto bbLeft =  m_context.bxLeft[faceIdx] * m_context.bxLeft[faceIdx] +
                       m_context.byLeft[faceIdx] * m_context.byLeft[faceIdx] +
                       m_context.bzLeft[faceIdx] * m_context.bzLeft[faceIdx];

        // Face-centered magnitude of the magnetic field on the right
        auto bbRight =  m_context.bxRight[faceIdx] * m_context.bxRight[faceIdx] +
                        m_context.byRight[faceIdx] * m_context.byRight[faceIdx] +
                        m_context.bzRight[faceIdx] * m_context.bzRight[faceIdx];
    
        // Face-centered normal velocity on the left
        auto uDotNLeft = m_context.uLeft[faceIdx] * m_context.faceNormalX[faceIdx] +
                         m_context.vLeft[faceIdx] * m_context.faceNormalY[faceIdx] +
                         m_context.wLeft[faceIdx] * m_context.faceNormalZ[faceIdx];

        // Face-centered normal velocity on the right
        auto uDotNRight = m_context.uRight[faceIdx] * m_context.faceNormalX[faceIdx] +
                          m_context.vRight[faceIdx] * m_context.faceNormalY[faceIdx] +
                          m_context.wRight[faceIdx] * m_context.faceNormalZ[faceIdx];

        // Face-centered normal magnetic field on the left
        auto bDotNLeft = m_context.bxLeft[faceIdx] * m_context.faceNormalX[faceIdx] +
                         m_context.byLeft[faceIdx] * m_context.faceNormalY[faceIdx] +
                         m_context.bzLeft[faceIdx] * m_context.faceNormalZ[faceIdx];

        // Face-centered normal magnetic field on the right
        auto bDotNRight = m_context.bxRight[faceIdx] * m_context.faceNormalX[faceIdx] +
                          m_context.byRight[faceIdx] * m_context.faceNormalY[faceIdx] +
                          m_context.bzRight[faceIdx] * m_context.faceNormalZ[faceIdx];
        
        // Face-centered dot product of magnetic field with velocity on the left
        auto bDotULeft = m_context.bxLeft[faceIdx] * m_context.uLeft[faceIdx] +
                         m_context.byLeft[faceIdx] * m_context.vLeft[faceIdx] +
                         m_context.bzLeft[faceIdx] * m_context.wLeft[faceIdx];
        
        // Face-centered dot product of magnetic field with velocity on the right
        auto bDotURight = m_context.bxRight[faceIdx] * m_context.uRight[faceIdx] +
                          m_context.byRight[faceIdx] * m_context.vRight[faceIdx] +
                          m_context.bzRight[faceIdx] * m_context.wRight[faceIdx];

        // Face-centered momentum densities on the left
        auto rhoULeft = m_context.rhoLeft[faceIdx] * m_context.uLeft[faceIdx];
        auto rhoVLeft = m_context.rhoLeft[faceIdx] * m_context.vLeft[faceIdx];
        auto rhoWLeft = m_context.rhoLeft[faceIdx] * m_context.wLeft[faceIdx];

        // Face-centered momentum densities on the right
        auto rhoURight = m_context.rhoRight[faceIdx] * m_context.uRight[faceIdx];
        auto rhoVRight = m_context.rhoRight[faceIdx] * m_context.vRight[faceIdx];
        auto rhoWRight = m_context.rhoRight[faceIdx] * m_context.wRight[faceIdx];

        // Face-centered total energy density on the left
        auto rhoELeft = m_context.rhoLeft[faceIdx] * (m_context.eLeft[faceIdx] + 0.5 * uuLeft);

        // Face-centered total energy density on the right
        auto rhoERight = m_context.rhoRight[faceIdx] * (m_context.eRight[faceIdx] + 0.5 * uuRight);

        // Local propagation speed
        double maxEigenValX = std::max(std::abs(m_context.uLeft[faceIdx]) + m_context.csLeft[faceIdx],
                                       std::abs(m_context.uRight[faceIdx]) + m_context.csRight[faceIdx]);
        double maxEigenValY = std::max(std::abs(m_context.vLeft[faceIdx]) + m_context.csLeft[faceIdx],
                                       std::abs(m_context.vRight[faceIdx]) + m_context.csRight[faceIdx]);
        double maxEigenValZ = std::max(std::abs(m_context.wLeft[faceIdx]) + m_context.csLeft[faceIdx],
                                       std::abs(m_context.wRight[faceIdx]) + m_context.csRight[faceIdx]);
        double maxEigenVal = maxEigenValX + maxEigenValY + maxEigenValZ;

        // Mass density flux
        m_context.rhoFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                               (m_context.rhoLeft[faceIdx] * uDotNLeft +
                                m_context.rhoRight[faceIdx] * uDotNRight -
                                maxEigenVal * (m_context.rhoRight[faceIdx] - m_context.rhoLeft[faceIdx]));

        // x-momentum density flux
        m_context.rhoUFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                (rhoULeft * uDotNLeft + m_context.pLeft[faceIdx] * m_context.faceNormalX[faceIdx] -
                                    m_context.bxLeft[faceIdx] * bDotNLeft +
                                rhoURight * uDotNRight + m_context.pRight[faceIdx] * m_context.faceNormalX[faceIdx] -
                                    m_context.bxRight[faceIdx] * bDotNRight -
                                maxEigenVal * (rhoURight - rhoULeft));
                                                             
        // y-momentum density flux
        m_context.rhoVFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                (rhoVRight * uDotNLeft + m_context.pLeft[faceIdx] * m_context.faceNormalY[faceIdx] -
                                    m_context.byLeft[faceIdx] * bDotNLeft +
                                rhoVRight * uDotNRight + m_context.pRight[faceIdx] * m_context.faceNormalY[faceIdx] -
                                    m_context.byRight[faceIdx] * bDotNRight -
                                maxEigenVal * (rhoVRight - rhoVLeft));

        // z-momentum density flux
        m_context.rhoWFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                (rhoWLeft * uDotNLeft + m_context.pLeft[faceIdx] * m_context.faceNormalZ[faceIdx] -
                                    m_context.bzLeft[faceIdx] * bDotNLeft +
                                rhoWRight * uDotNRight + m_context.pRight[faceIdx] * m_context.faceNormalZ[faceIdx] -
                                    m_context.bzRight[faceIdx] * bDotNRight -
                                maxEigenVal * (rhoWRight - rhoWLeft));

        // total energy density flux
        m_context.rhoEFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                ((rhoELeft + m_context.pLeft[faceIdx]) * uDotNLeft - bDotULeft * bDotNLeft +
                                (rhoERight + m_context.pRight[faceIdx]) * uDotNRight -bDotURight * bDotNRight -
                                maxEigenVal * (rhoERight - rhoELeft));

        // magnetic field fluxes
        m_context.bxFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                    ((m_context.uLeft[faceIdx] * bDotNLeft - m_context.bxLeft[faceIdx] * uDotNLeft) -
                                     (m_context.uRight[faceIdx] * bDotNRight - m_context.bxRight[faceIdx] * uDotNRight) -
                                     maxEigenVal * (m_context.bxRight[faceIdx] - m_context.bxLeft[faceIdx]));
        m_context.byFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                     ((m_context.vLeft[faceIdx] * bDotNLeft - m_context.byLeft[faceIdx] * uDotNLeft) -
                                      (m_context.vRight[faceIdx] * bDotNRight - m_context.byRight[faceIdx] * uDotNRight) -
                                      maxEigenVal * (m_context.byRight[faceIdx] - m_context.byLeft[faceIdx]));

        m_context.bzFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                      ((m_context.wLeft[faceIdx] * bDotNLeft - m_context.bzLeft[faceIdx] * uDotNLeft) -
                                       (m_context.wRight[faceIdx] * bDotNRight - m_context.bzRight[faceIdx] * uDotNRight) -
                                       maxEigenVal * (m_context.bzRight[faceIdx] - m_context.bzLeft[faceIdx]));
    }

    FluxContext& m_context;
};

struct GodunovConstantFluxKernel {
    GodunovConstantFluxKernel(FluxContext& context) : m_context(context) {}

    inline void operator()(std::size_t const i) {
        std::size_t const faceIdx = m_context.faceIdxs[i];

        // Face-centered magnitude of the velocity on the left
        auto uuLeft =  m_context.uLeft[faceIdx] * m_context.uLeft[faceIdx] +
                       m_context.vLeft[faceIdx] * m_context.vLeft[faceIdx] +
                       m_context.wLeft[faceIdx] * m_context.wLeft[faceIdx];
    
        // Face-centered normal velocity on the left
        auto uDotNLeft = m_context.uLeft[faceIdx] * m_context.faceNormalX[faceIdx] +
                         m_context.vLeft[faceIdx] * m_context.faceNormalY[faceIdx] +
                         m_context.wLeft[faceIdx] * m_context.faceNormalZ[faceIdx];

        // Face-centered momentum densities on the left
        auto rhoULeft = m_context.rhoLeft[faceIdx] * m_context.uLeft[faceIdx];
        auto rhoVLeft = m_context.rhoLeft[faceIdx] * m_context.vLeft[faceIdx];
        auto rhoWLeft = m_context.rhoLeft[faceIdx] * m_context.wLeft[faceIdx];

        // Face-centered total energy density on the left
        auto rhoELeft = m_context.rhoLeft[faceIdx] * (m_context.eLeft[faceIdx] + 0.5 * uuLeft);

        // Mass density flux
        m_context.rhoFlux[faceIdx] = m_context.faceArea[faceIdx] *
                                     (m_context.rhoLeft[faceIdx] * uDotNLeft);

        // x-momentum density flux
        m_context.rhoUFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                      (rhoULeft * uDotNLeft + m_context.pLeft[faceIdx] * m_context.faceNormalX[faceIdx]);
                                                             
        // y-momentum density flux
        m_context.rhoVFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                      (rhoVLeft * uDotNLeft + m_context.pLeft[faceIdx] * m_context.faceNormalY[faceIdx]);

        // z-momentum density flux
        m_context.rhoWFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                      (rhoWLeft * uDotNLeft + m_context.pLeft[faceIdx] * m_context.faceNormalZ[faceIdx]);

        // total energy density flux
        m_context.rhoEFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                      ((rhoELeft + m_context.pLeft[faceIdx]) * uDotNLeft);
    }

    FluxContext& m_context;
};

FluxContext::FluxContext(IGrid const& grid, ReconstructionContext const& rc) :
    numFaces(grid.NumFaces()), faceIdxToNodeIdxs(grid.FaceIdxToCellIdxs()), faceArea(grid.FaceAreas()),
    faceNormalX(grid.FaceNormalX()), faceNormalY(grid.FaceNormalY()), faceNormalZ(grid.FaceNormalZ()), faceIdxs(grid.FaceIdxs()),
    rhoLeft(rc.rhoLeft), uLeft(rc.uLeft), vLeft(rc.vLeft), wLeft(rc.wLeft), pLeft(rc.pLeft), eLeft(rc.eLeft), csLeft(rc.csLeft),
    rhoRight(rc.rhoRight), uRight(rc.uRight), vRight(rc.vRight), wRight(rc.wRight), pRight(rc.pRight), eRight(rc.eRight), csRight(rc.csRight),
    bxLeft(rc.bxLeft), byLeft(rc.byLeft), bzLeft(rc.bzLeft), bxRight(rc.bxRight), byRight(rc.byRight), bzRight(rc.bzRight) {
    rhoFlux.resize(numFaces, 0.0);
    rhoUFlux.resize(numFaces, 0.0);
    rhoVFlux.resize(numFaces, 0.0);
    rhoWFlux.resize(numFaces, 0.0);
    rhoEFlux.resize(numFaces, 0.0);
    bxFlux.resize(numFaces, 0.0);
    byFlux.resize(numFaces, 0.0);
    bzFlux.resize(numFaces, 0.0);
}

class GodunovConstantFlux : public IFlux {
public:
    GodunovConstantFlux(IGrid const& grid, ReconstructionContext const& rc) {
        m_context = std::make_unique<FluxContext>(grid, rc);
    }

    void ComputeInterfaceFluxes(ExecutionController const& execCtrl) const {
        GodunovConstantFluxKernel kern(*m_context);
        execCtrl.LaunchKernel(kern, m_context->numFaces);
    }
};

class KTFlux : public IFlux {
public:
    KTFlux(IGrid const& grid, ReconstructionContext const& rc) {
        m_context = std::make_unique<FluxContext>(grid, rc);
    }

    void ComputeInterfaceFluxes(ExecutionController const& execCtrl) const {
        KTFluxKernel kern(*m_context);
        execCtrl.LaunchKernel(kern, m_context->numFaces);
    }
};

std::unique_ptr<IFlux> fluxFactory(Profile const& profile, IGrid const& grid, ReconstructionContext const& rc) {
    if (FluxScheme::KT == profile.m_fluxOption) {
        return std::make_unique<KTFlux>(grid, rc);
    }
    throw Error::INVALID_FLUX_SCHEME;
}

} // namespace MHD
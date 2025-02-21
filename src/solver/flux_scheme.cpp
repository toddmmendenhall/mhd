#include <context.hpp>
#include <constants.hpp>
#include <error.hpp>
#include <execution_controller.hpp>
#include <flux_scheme.hpp>
#include <profile.hpp>
#include <profile_options.hpp>
#include <face.hpp>

namespace MHD {

struct HLLCKernel {
    HLLCKernel(FluxContext& fluxContext) : fluxContext(fluxContext) {}

    inline void operator()(std::size_t const idx) {}

    FluxContext& fluxContext;
};

struct KTKernel {
    KTKernel(FluxContext& fluxContext) : fluxContext(fluxContext) {}

    inline void operator()(std::size_t const idx) {}

    FluxContext& fluxContext;
};

struct HighOrderGodunovKernel {
    HighOrderGodunovKernel(FluxContext& fluxContext) : m_context(fluxContext) {}

    inline void operator()(std::size_t const faceIdx) {
        // normal component of the velocity on the left of the face
        auto uBarLeft = m_context.uLeft[faceIdx] * m_context.faceNormalX[faceIdx] +
                        m_context.vLeft[faceIdx] * m_context.faceNormalY[faceIdx] +
                        m_context.wLeft[faceIdx] * m_context.faceNormalZ[faceIdx];

        // normal component of the velocity on the right of the face
        auto uBarRight = m_context.uRight[faceIdx] * m_context.faceNormalX[faceIdx] +
                         m_context.vRight[faceIdx] * m_context.faceNormalY[faceIdx] +
                         m_context.wRight[faceIdx] * m_context.faceNormalZ[faceIdx];

        // normal component of the magnetic field on the left of the face
        auto bBarLeft = m_context.bXLeft[faceIdx] * m_context.faceNormalX[faceIdx] +
                        m_context.bYLeft[faceIdx] * m_context.faceNormalY[faceIdx] +
                        m_context.bZLeft[faceIdx] * m_context.faceNormalZ[faceIdx];

        // normal component of the magnetic field on the right of the face
        auto bBarRight = m_context.bXRight[faceIdx] * m_context.faceNormalX[faceIdx] +
                         m_context.bYRight[faceIdx] * m_context.faceNormalY[faceIdx] +
                         m_context.bZRight[faceIdx] * m_context.faceNormalZ[faceIdx];
        
        // cell-centered momentum densities on the left of the face
        auto rhoULeft = m_context.rhoLeft[faceIdx] * m_context.uLeft[faceIdx];
        auto rhoVLeft = m_context.rhoLeft[faceIdx] * m_context.vLeft[faceIdx];
        auto rhoWLeft = m_context.rhoLeft[faceIdx] * m_context.wLeft[faceIdx];
        
        // cell-centered momentum densities on the right of the face
        auto rhoURight = m_context.rhoRight[faceIdx] * m_context.uRight[faceIdx];
        auto rhoVRight = m_context.rhoRight[faceIdx] * m_context.vRight[faceIdx];
        auto rhoWRight = m_context.rhoRight[faceIdx] * m_context.wRight[faceIdx];

        // cell-centered total energy density on the left of the face
        auto rhoELeft = m_context.rhoLeft[faceIdx] * (m_context.eLeft[faceIdx] + 0.5 * m_context.uuLeft[faceIdx]) + 0.5 * VACUUM_PERMEABILITY_INV * m_context.bbLeft[faceIdx];

        // cell-centered total energy density on the right of the face
        auto rhoERight = m_context.rhoRight[faceIdx] * (m_context.eRight[faceIdx] + 0.5 * m_context.uuRight[faceIdx]) + 0.5 * VACUUM_PERMEABILITY_INV * m_context.bbRight[faceIdx];

        // cell-centered dot(u,B) on the left of the face
        auto uDotBLeft = m_context.uLeft[faceIdx] * m_context.bXLeft[faceIdx] +
                         m_context.vLeft[faceIdx] * m_context.bYLeft[faceIdx] +
                         m_context.wLeft[faceIdx] * m_context.bZLeft[faceIdx];
        
        // cell-centered dot(u,B) on the right of the face
        auto uDotBRight = m_context.uRight[faceIdx] * m_context.bXRight[faceIdx] +
                          m_context.vRight[faceIdx] * m_context.bYRight[faceIdx] +
                          m_context.wRight[faceIdx] * m_context.bZRight[faceIdx];
        
        // face-centered mass density flux
        m_context.rhoFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                     (m_context.rhoLeft[faceIdx] * uBarLeft +
                                     m_context.rhoRight[faceIdx] * uBarRight);

        // face-centered x-momentum density flux
        m_context.rhoUFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                      (rhoULeft * uBarLeft + m_context.pLeft[faceIdx] * m_context.faceNormalX[faceIdx] - m_context.bXLeft[faceIdx] * bBarLeft +
                                      rhoURight * uBarRight + m_context.pRight[faceIdx] * m_context.faceNormalX[faceIdx] - m_context.bXRight[faceIdx] * bBarRight);
                                                             
        // face-centered y-momentum density flux
        m_context.rhoVFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                      (rhoVLeft * uBarLeft + m_context.pLeft[faceIdx] * m_context.faceNormalY[faceIdx] - m_context.bYLeft[faceIdx] * bBarLeft +
                                      rhoVRight * uBarRight + m_context.pRight[faceIdx] * m_context.faceNormalY[faceIdx] - m_context.bYRight[faceIdx] * bBarRight);
        
        // face-centered z-momentum density flux
        m_context.rhoWFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                      (rhoWLeft * uBarLeft + m_context.pLeft[faceIdx] * m_context.faceNormalZ[faceIdx] - m_context.bZLeft[faceIdx] * bBarLeft +
                                      rhoWRight * uBarRight + m_context.pRight[faceIdx] * m_context.faceNormalZ[faceIdx] - m_context.bZRight[faceIdx] * bBarRight);
        
        // face-centered total energy density flux
        m_context.rhoEFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                      ((rhoELeft + m_context.pLeft[faceIdx]) * uBarLeft - uDotBLeft * bBarLeft +
                                      (rhoERight + m_context.pRight[faceIdx]) * uBarRight - uDotBRight * bBarRight);
        
        // face-centered x-magnetic flux
        m_context.bXFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                    (m_context.uLeft[faceIdx] * bBarLeft - m_context.bXLeft[faceIdx] * uBarLeft +
                                    m_context.uRight[faceIdx] * bBarRight - m_context.bXRight[faceIdx] * uBarRight);

        // face-centered y-magnetic flux
        m_context.bYFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                    (m_context.vLeft[faceIdx] * bBarLeft - m_context.bYLeft[faceIdx] * uBarLeft +
                                    m_context.vRight[faceIdx] * bBarRight - m_context.bYRight[faceIdx] * uBarRight);

        // face-centered z-magnetic flux
        m_context.bZFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                    (m_context.wLeft[faceIdx] * bBarLeft - m_context.bZLeft[faceIdx] * uBarLeft +
                                    m_context.wRight[faceIdx] * bBarRight - m_context.bZRight[faceIdx] * uBarRight);
    }

    FluxContext& m_context;
};

struct LowOrderGodunovKernel {
    LowOrderGodunovKernel(FluxContext& context) : m_context(context) {}

    inline void operator()(std::size_t const i) {
        // normal component of the velocity inside the cell
        auto uBar = m_context.uLeft[i] * m_context.faceNormalX[i] +
                    m_context.vLeft[i] * m_context.faceNormalY[i] +
                    m_context.wLeft[i] * m_context.faceNormalZ[i];

        // normal component of the magnetic field inside the cell
        auto bBar = m_context.bXLeft[i] * m_context.faceNormalX[i] +
                    m_context.bYLeft[i] * m_context.faceNormalY[i] +
                    m_context.bZLeft[i] * m_context.faceNormalZ[i];
        
        // cell-centered momentum densities inside the cell
        auto rhoU = m_context.rhoLeft[i] * m_context.uLeft[i];
        auto rhoV = m_context.rhoLeft[i] * m_context.vLeft[i];
        auto rhoW = m_context.rhoLeft[i] * m_context.wLeft[i];

        // cell-centered total energy density
        auto rhoE = m_context.rhoLeft[i] * (m_context.eLeft[i] + 0.5 * m_context.uuLeft[i]) + 0.5 * VACUUM_PERMEABILITY_INV * m_context.bbLeft[i];

        // cell-centered dot(u,B)
        auto uDotB = m_context.uLeft[i] * m_context.bXLeft[i] +
                     m_context.vLeft[i] * m_context.bYLeft[i] +
                     m_context.wLeft[i] * m_context.bZLeft[i];

        // cell-centered mass density flux
        m_context.rhoFlux[i] = m_context.faceArea[i] * m_context.rhoLeft[i] * uBar;

        // cell-centered x-momentum density flux
        m_context.rhoUFlux[i] = m_context.faceArea[i] * (rhoU * uBar + m_context.pLeft[i] * m_context.faceNormalX[i] - m_context.bXLeft[i] * bBar);
                                                             
        // cell-centered y-momentum density flux
        m_context.rhoVFlux[i] = m_context.faceArea[i] * (rhoV * uBar + m_context.pLeft[i] * m_context.faceNormalY[i] - m_context.bYLeft[i] * bBar);

        // cell-centered z-momentum density flux
        m_context.rhoWFlux[i] = m_context.faceArea[i] * (rhoW * uBar + m_context.pLeft[i] * m_context.faceNormalZ[i] - m_context.bZLeft[i] * bBar);

        // cell-centered total energy density flux
        m_context.rhoEFlux[i] = m_context.faceArea[i] * ((rhoE + m_context.pLeft[i]) * uBar - uDotB * bBar);
        
        // cell-centered x-magnetic flux
        m_context.bXFlux[i] = m_context.faceArea[i] * (m_context.uLeft[i] * bBar - m_context.bXLeft[i] * uBar);

        // cell-centered y-magnetic flux
        m_context.bYFlux[i] = m_context.faceArea[i] * (m_context.vLeft[i] * bBar - m_context.bYLeft[i] * uBar);

        // cell-centered z-magnetic flux
        m_context.bZFlux[i] = m_context.faceArea[i] * (m_context.wLeft[i] * bBar - m_context.bZLeft[i] * uBar);
    }

    FluxContext& m_context;
};

class HLLC : public IFluxScheme {
public:
    void ComputeInterfaceFluxes(ExecutionController const& execCtrl, FluxContext& fluxContext) const {
        HLLCKernel kern(fluxContext);
        execCtrl.LaunchKernel(kern, fluxContext.faceArea.size());
    }
};

class KT : public IFluxScheme {
public:
    void ComputeInterfaceFluxes(ExecutionController const& execCtrl, FluxContext& fluxContext) const {
        KTKernel kern(fluxContext);
        execCtrl.LaunchKernel(kern, fluxContext.faceArea.size());
    }
};

class HighOrderGodunov : public IFluxScheme {
public:
    void ComputeInterfaceFluxes(ExecutionController const& execCtrl, FluxContext& fluxContext) const {
        HighOrderGodunovKernel kern(fluxContext);
        execCtrl.LaunchKernel(kern, fluxContext.faceArea.size());
    }
};

class LowOrderGodunov : public IFluxScheme {
public:
    void ComputeInterfaceFluxes(ExecutionController const& execCtrl, FluxContext& fluxContext) const {
        LowOrderGodunovKernel kern(fluxContext);
        execCtrl.LaunchKernel(kern, fluxContext.faceArea.size());
    }
};

std::unique_ptr<IFluxScheme> fluxSchemeFactory(Profile const& profile) {
    if (FluxScheme::HLLC == profile.m_fluxOption) {
        return std::make_unique<HLLC>();
    }
    if (FluxScheme::KT == profile.m_fluxOption) {
        return std::make_unique<KT>();
    }
    if (FluxScheme::HOG == profile.m_fluxOption) {
        return std::make_unique<HighOrderGodunov>();
    }
    if (FluxScheme::LOG == profile.m_fluxOption) {
        return std::make_unique<LowOrderGodunov>();
    }
    throw Error::INVALID_FLUX_SCHEME;
}

} // namespace MHD
#include <context.hpp>
#include <constants.hpp>
#include <error.hpp>
#include <execution_controller.hpp>
#include <flux_scheme.hpp>
#include <profile.hpp>
#include <profile_options.hpp>
#include <face.hpp>

namespace MHD {

// struct HighOrderGodunovKernel {
//     HighOrderGodunovKernel(FluxContext& fluxContext) : m_context(fluxContext) {}

//     inline void operator()(std::size_t const faceIdx) {
//         // normal component of the velocity on the left of the face
//         auto uBarLeft = m_context.uLeft[faceIdx] * m_context.faceNormalX[faceIdx] +
//                         m_context.vLeft[faceIdx] * m_context.faceNormalY[faceIdx] +
//                         m_context.wLeft[faceIdx] * m_context.faceNormalZ[faceIdx];

//         // normal component of the velocity on the right of the face
//         auto uBarRight = m_context.uRight[faceIdx] * m_context.faceNormalX[faceIdx] +
//                          m_context.vRight[faceIdx] * m_context.faceNormalY[faceIdx] +
//                          m_context.wRight[faceIdx] * m_context.faceNormalZ[faceIdx];

//         // normal component of the magnetic field on the left of the face
//         auto bBarLeft = m_context.bXLeft[faceIdx] * m_context.faceNormalX[faceIdx] +
//                         m_context.bYLeft[faceIdx] * m_context.faceNormalY[faceIdx] +
//                         m_context.bZLeft[faceIdx] * m_context.faceNormalZ[faceIdx];

//         // normal component of the magnetic field on the right of the face
//         auto bBarRight = m_context.bXRight[faceIdx] * m_context.faceNormalX[faceIdx] +
//                          m_context.bYRight[faceIdx] * m_context.faceNormalY[faceIdx] +
//                          m_context.bZRight[faceIdx] * m_context.faceNormalZ[faceIdx];
        
//         // cell-centered momentum densities on the left of the face
//         auto rhoULeft = m_context.rhoLeft[faceIdx] * m_context.uLeft[faceIdx];
//         auto rhoVLeft = m_context.rhoLeft[faceIdx] * m_context.vLeft[faceIdx];
//         auto rhoWLeft = m_context.rhoLeft[faceIdx] * m_context.wLeft[faceIdx];
        
//         // cell-centered momentum densities on the right of the face
//         auto rhoURight = m_context.rhoRight[faceIdx] * m_context.uRight[faceIdx];
//         auto rhoVRight = m_context.rhoRight[faceIdx] * m_context.vRight[faceIdx];
//         auto rhoWRight = m_context.rhoRight[faceIdx] * m_context.wRight[faceIdx];

//         // cell-centered total energy density on the left of the face
//         auto rhoELeft = m_context.rhoLeft[faceIdx] * (m_context.eLeft[faceIdx] + 0.5 * m_context.uuLeft[faceIdx]) + 0.5 * VACUUM_PERMEABILITY_INV * m_context.bbLeft[faceIdx];

//         // cell-centered total energy density on the right of the face
//         auto rhoERight = m_context.rhoRight[faceIdx] * (m_context.eRight[faceIdx] + 0.5 * m_context.uuRight[faceIdx]) + 0.5 * VACUUM_PERMEABILITY_INV * m_context.bbRight[faceIdx];

//         // cell-centered dot(u,B) on the left of the face
//         auto uDotBLeft = m_context.uLeft[faceIdx] * m_context.bXLeft[faceIdx] +
//                          m_context.vLeft[faceIdx] * m_context.bYLeft[faceIdx] +
//                          m_context.wLeft[faceIdx] * m_context.bZLeft[faceIdx];
        
//         // cell-centered dot(u,B) on the right of the face
//         auto uDotBRight = m_context.uRight[faceIdx] * m_context.bXRight[faceIdx] +
//                           m_context.vRight[faceIdx] * m_context.bYRight[faceIdx] +
//                           m_context.wRight[faceIdx] * m_context.bZRight[faceIdx];
        
//         // face-centered mass density flux
//         m_context.rhoFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
//                                      (m_context.rhoLeft[faceIdx] * uBarLeft +
//                                      m_context.rhoRight[faceIdx] * uBarRight);

//         // face-centered x-momentum density flux
//         m_context.rhoUFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
//                                       (rhoULeft * uBarLeft + m_context.pLeft[faceIdx] * m_context.faceNormalX[faceIdx] - m_context.bXLeft[faceIdx] * bBarLeft +
//                                       rhoURight * uBarRight + m_context.pRight[faceIdx] * m_context.faceNormalX[faceIdx] - m_context.bXRight[faceIdx] * bBarRight);
                                                             
//         // face-centered y-momentum density flux
//         m_context.rhoVFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
//                                       (rhoVLeft * uBarLeft + m_context.pLeft[faceIdx] * m_context.faceNormalY[faceIdx] - m_context.bYLeft[faceIdx] * bBarLeft +
//                                       rhoVRight * uBarRight + m_context.pRight[faceIdx] * m_context.faceNormalY[faceIdx] - m_context.bYRight[faceIdx] * bBarRight);
        
//         // face-centered z-momentum density flux
//         m_context.rhoWFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
//                                       (rhoWLeft * uBarLeft + m_context.pLeft[faceIdx] * m_context.faceNormalZ[faceIdx] - m_context.bZLeft[faceIdx] * bBarLeft +
//                                       rhoWRight * uBarRight + m_context.pRight[faceIdx] * m_context.faceNormalZ[faceIdx] - m_context.bZRight[faceIdx] * bBarRight);
        
//         // face-centered total energy density flux
//         m_context.rhoEFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
//                                       ((rhoELeft + m_context.pLeft[faceIdx]) * uBarLeft - uDotBLeft * bBarLeft +
//                                       (rhoERight + m_context.pRight[faceIdx]) * uBarRight - uDotBRight * bBarRight);
        
//         // face-centered x-magnetic flux
//         m_context.bXFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
//                                     (m_context.uLeft[faceIdx] * bBarLeft - m_context.bXLeft[faceIdx] * uBarLeft +
//                                     m_context.uRight[faceIdx] * bBarRight - m_context.bXRight[faceIdx] * uBarRight);

//         // face-centered y-magnetic flux
//         m_context.bYFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
//                                     (m_context.vLeft[faceIdx] * bBarLeft - m_context.bYLeft[faceIdx] * uBarLeft +
//                                     m_context.vRight[faceIdx] * bBarRight - m_context.bYRight[faceIdx] * uBarRight);

//         // face-centered z-magnetic flux
//         m_context.bZFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
//                                     (m_context.wLeft[faceIdx] * bBarLeft - m_context.bZLeft[faceIdx] * uBarLeft +
//                                     m_context.wRight[faceIdx] * bBarRight - m_context.bZRight[faceIdx] * uBarRight);
//     }

//     FluxContext& m_context;
// };

struct LowOrderGodunovKernel {
    LowOrderGodunovKernel(FluxContext& context) : m_context(context) {}

    inline void operator()(std::size_t const i) {
        // Face-centered magnitude of the velocity 
        auto uu =  m_context.uFace[i] * m_context.uFace[i] +
                   m_context.vFace[i] * m_context.vFace[i] +
                   m_context.wFace[i] * m_context.wFace[i];
    
        // Face-centered normal velocity
        auto uBar = m_context.uFace[i] * m_context.faceNormalX[i] +
                    m_context.vFace[i] * m_context.faceNormalY[i] +
                    m_context.wFace[i] * m_context.faceNormalZ[i];

        // Face-centered momentum densities
        auto rhoU = m_context.rhoFace[i] * m_context.uFace[i];
        auto rhoV = m_context.rhoFace[i] * m_context.vFace[i];
        auto rhoW = m_context.rhoFace[i] * m_context.wFace[i];

        // Face-centered total energy density
        auto rhoE = m_context.rhoFace[i] * (m_context.eFace[i] + 0.5 * uu);

        // Mass density flux
        m_context.rhoFlux[i] = m_context.faceArea[i] * m_context.rhoFace[i] * uBar;

        // x-momentum density flux
        m_context.rhoUFlux[i] = m_context.faceArea[i] * (rhoU * uBar + m_context.pFace[i] * m_context.faceNormalX[i]);
                                                             
        // y-momentum density flux
        m_context.rhoVFlux[i] = m_context.faceArea[i] * (rhoV * uBar + m_context.pFace[i] * m_context.faceNormalY[i]);

        // z-momentum density flux
        m_context.rhoWFlux[i] = m_context.faceArea[i] * (rhoW * uBar + m_context.pFace[i] * m_context.faceNormalZ[i]);

        // total energy density flux
        m_context.rhoEFlux[i] = m_context.faceArea[i] * ((rhoE + m_context.pFace[i]) * uBar);
    }

    FluxContext& m_context;
};

class LowOrderGodunov : public IFluxScheme {
public:
    LowOrderGodunov(FluxContext& context) :
        m_context(context) {}

    void ComputeInterfaceFluxes(ExecutionController const& execCtrl) const {
        LowOrderGodunovKernel kern(m_context);
        execCtrl.LaunchKernel(kern, m_context.numFaces);
    }

    FluxContext& m_context;
};

std::unique_ptr<IFluxScheme> fluxSchemeFactory(Profile const& profile, FluxContext& context) {
    if (FluxScheme::LOG == profile.m_fluxOption) {
        return std::make_unique<LowOrderGodunov>(context);
    }
    throw Error::INVALID_FLUX_SCHEME;
}

} // namespace MHD
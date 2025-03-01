#include <context.hpp>
#include <constants.hpp>
#include <error.hpp>
#include <execution_controller.hpp>
#include <flux/flux_scheme.hpp>
#include <profile.hpp>
#include <profile_options.hpp>

#include <cstdlib>

namespace MHD {

struct LowOrderGodunovKernel {
    LowOrderGodunovKernel(FluxContext& context) : m_context(context) {}

    inline void operator()(std::size_t const faceIdx) {
        // Face-centered magnitude of the velocity on the left
        auto uuLeft =  m_context.uLeft[faceIdx] * m_context.uLeft[faceIdx] +
                       m_context.vLeft[faceIdx] * m_context.vLeft[faceIdx] +
                       m_context.wLeft[faceIdx] * m_context.wLeft[faceIdx];

        // Face-centered magnitude of the velocity on the right
        auto uuRight =  m_context.uRight[faceIdx] * m_context.uRight[faceIdx] +
                        m_context.vRight[faceIdx] * m_context.vRight[faceIdx] +
                        m_context.wRight[faceIdx] * m_context.wRight[faceIdx];
    
        // Face-centered normal velocity on the left
        auto uDotNLeft = m_context.uLeft[faceIdx] * m_context.faceNormalX[faceIdx] +
                         m_context.vLeft[faceIdx] * m_context.faceNormalY[faceIdx] +
                         m_context.wLeft[faceIdx] * m_context.faceNormalZ[faceIdx];

        // Face-centered normal velocity on the right
        auto uDotNRight = m_context.uRight[faceIdx] * m_context.faceNormalX[faceIdx] +
                          m_context.vRight[faceIdx] * m_context.faceNormalY[faceIdx] +
                          m_context.wRight[faceIdx] * m_context.faceNormalZ[faceIdx];

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
        double maxEigenValX = std::max(std::max(std::abs(m_context.uLeft[faceIdx] - m_context.csLeft[faceIdx]), std::abs(m_context.uLeft[faceIdx] + m_context.csLeft[faceIdx])),
                                       std::max(std::abs(m_context.uRight[faceIdx] - m_context.csRight[faceIdx]), std::abs(m_context.uRight[faceIdx] + m_context.csRight[faceIdx])));
        double maxEigenValY = std::max(std::max(std::abs(m_context.vLeft[faceIdx] - m_context.csLeft[faceIdx]), std::abs(m_context.vLeft[faceIdx] + m_context.csLeft[faceIdx])),
                                       std::max(std::abs(m_context.vRight[faceIdx] - m_context.csRight[faceIdx]), std::abs(m_context.vRight[faceIdx] + m_context.csRight[faceIdx])));
        double maxEigenValZ = std::max(std::max(std::abs(m_context.wLeft[faceIdx] - m_context.csLeft[faceIdx]), std::abs(m_context.wLeft[faceIdx] + m_context.csLeft[faceIdx])),
                                       std::max(std::abs(m_context.wRight[faceIdx] - m_context.csRight[faceIdx]), std::abs(m_context.wRight[faceIdx] + m_context.csRight[faceIdx])));
        double maxEigenVal = maxEigenValX + maxEigenValY + maxEigenValZ;

        // Mass density flux
        m_context.rhoFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                               (m_context.rhoLeft[faceIdx] * uDotNLeft +
                                m_context.rhoRight[faceIdx] * uDotNRight -
                                maxEigenVal * (m_context.rhoRight[faceIdx] - m_context.rhoLeft[faceIdx]));

        // x-momentum density flux
        m_context.rhoUFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                (rhoULeft * uDotNLeft + m_context.pLeft[faceIdx] * m_context.faceNormalX[faceIdx] +
                                rhoURight * uDotNRight + m_context.pRight[faceIdx] * -m_context.faceNormalX[faceIdx] -
                                maxEigenVal * (rhoURight - rhoULeft));
                                                             
        // y-momentum density flux
        m_context.rhoVFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                (rhoVRight * uDotNLeft + m_context.pLeft[faceIdx] * m_context.faceNormalY[faceIdx] +
                                rhoVRight * uDotNRight + m_context.pRight[faceIdx] * -m_context.faceNormalY[faceIdx] -
                                maxEigenVal * (rhoVRight - rhoVLeft));

        // z-momentum density flux
        m_context.rhoWFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                (rhoWLeft * uDotNLeft + m_context.pLeft[faceIdx] * m_context.faceNormalZ[faceIdx] +
                                rhoWRight * uDotNRight + m_context.pRight[faceIdx] * -m_context.faceNormalZ[faceIdx] -
                                maxEigenVal * (rhoWRight - rhoWLeft));

        // total energy density flux
        m_context.rhoEFlux[faceIdx] = 0.5 * m_context.faceArea[faceIdx] *
                                ((rhoELeft + m_context.pLeft[faceIdx]) * uDotNLeft +
                                (rhoERight + m_context.pRight[faceIdx]) * uDotNRight -
                                maxEigenVal * (rhoERight - rhoELeft));
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
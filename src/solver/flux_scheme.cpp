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
    HighOrderGodunovKernel(FluxContext& fluxContext) : ctx(fluxContext) {}

    inline void operator()(std::size_t const faceIdx) {
        // normal component of the velocity on the left of the face
        auto uBarLeft = ctx.uLeft[faceIdx] * ctx.faceNormalX[faceIdx] +
                        ctx.vLeft[faceIdx] * ctx.faceNormalY[faceIdx] +
                        ctx.wLeft[faceIdx] * ctx.faceNormalZ[faceIdx];

        // normal component of the velocity on the right of the face
        auto uBarRight = ctx.uRight[faceIdx] * ctx.faceNormalX[faceIdx] +
                         ctx.vRight[faceIdx] * ctx.faceNormalY[faceIdx] +
                         ctx.wRight[faceIdx] * ctx.faceNormalZ[faceIdx];

        // normal component of the magnetic field on the left of the face
        auto bBarLeft = ctx.bXLeft[faceIdx] * ctx.faceNormalX[faceIdx] +
                        ctx.bYLeft[faceIdx] * ctx.faceNormalY[faceIdx] +
                        ctx.bZLeft[faceIdx] * ctx.faceNormalZ[faceIdx];

        // normal component of the magnetic field on the right of the face
        auto bBarRight = ctx.bXRight[faceIdx] * ctx.faceNormalX[faceIdx] +
                         ctx.bYRight[faceIdx] * ctx.faceNormalY[faceIdx] +
                         ctx.bZRight[faceIdx] * ctx.faceNormalZ[faceIdx];
        
        // cell-centered momentum densities on the left of the face
        auto rhoULeft = ctx.rhoLeft[faceIdx] * ctx.uLeft[faceIdx];
        auto rhoVLeft = ctx.rhoLeft[faceIdx] * ctx.vLeft[faceIdx];
        auto rhoWLeft = ctx.rhoLeft[faceIdx] * ctx.wLeft[faceIdx];
        
        // cell-centered momentum densities on the right of the face
        auto rhoURight = ctx.rhoRight[faceIdx] * ctx.uRight[faceIdx];
        auto rhoVRight = ctx.rhoRight[faceIdx] * ctx.vRight[faceIdx];
        auto rhoWRight = ctx.rhoRight[faceIdx] * ctx.wRight[faceIdx];

        // cell-centered total energy density on the left of the face
        auto rhoELeft = ctx.rhoLeft[faceIdx] * (ctx.eLeft[faceIdx] + 0.5 * ctx.uuLeft[faceIdx]) + 0.5 * VACUUM_PERMITTIVITY_INV * ctx.bbLeft[faceIdx];

        // cell-centered total energy density on the right of the face
        auto rhoERight = ctx.rhoRight[faceIdx] * (ctx.eRight[faceIdx] + 0.5 * ctx.uuRight[faceIdx]) + 0.5 * VACUUM_PERMITTIVITY_INV * ctx.bbRight[faceIdx];

        // cell-centered dot(u,B) on the left of the face
        auto uDotBLeft = ctx.uLeft[faceIdx] * ctx.bXLeft[faceIdx] +
                         ctx.vLeft[faceIdx] * ctx.bYLeft[faceIdx] +
                         ctx.wLeft[faceIdx] * ctx.bZLeft[faceIdx];
        
        // cell-centered dot(u,B) on the right of the face
        auto uDotBRight = ctx.uRight[faceIdx] * ctx.bXRight[faceIdx] +
                          ctx.vRight[faceIdx] * ctx.bYRight[faceIdx] +
                          ctx.wRight[faceIdx] * ctx.bZRight[faceIdx];
        
        // face-centered mass density flux
        ctx.rhoFlux[faceIdx] = 0.5 * ctx.faceArea[faceIdx] *
                               (ctx.rhoLeft[faceIdx] * uBarLeft +
                               ctx.rhoRight[faceIdx] * uBarRight);

        // face-centered x-momentum density flux
        ctx.rhoUFlux[faceIdx] = 0.5 * ctx.faceArea[faceIdx] *
                                (rhoULeft * uBarLeft + ctx.pLeft[faceIdx] * ctx.faceNormalX[faceIdx] - ctx.bXLeft[faceIdx] * bBarLeft +
                                rhoURight * uBarRight + ctx.pRight[faceIdx] * ctx.faceNormalX[faceIdx] - ctx.bXRight[faceIdx] * bBarRight);
                                                             
        // face-centered y-momentum density flux
        ctx.rhoVFlux[faceIdx] = 0.5 * ctx.faceArea[faceIdx] *
                                (rhoVLeft * uBarLeft + ctx.pLeft[faceIdx] * ctx.faceNormalY[faceIdx] - ctx.bYLeft[faceIdx] * bBarLeft +
                                rhoVRight * uBarRight + ctx.pRight[faceIdx] * ctx.faceNormalY[faceIdx] - ctx.bYRight[faceIdx] * bBarRight);
        
        // face-centered z-momentum density flux
        ctx.rhoWFlux[faceIdx] = 0.5 * ctx.faceArea[faceIdx] *
                                (rhoWLeft * uBarLeft + ctx.pLeft[faceIdx] * ctx.faceNormalZ[faceIdx] - ctx.bZLeft[faceIdx] * bBarLeft +
                                rhoWRight * uBarRight + ctx.pRight[faceIdx] * ctx.faceNormalZ[faceIdx] - ctx.bZRight[faceIdx] * bBarRight);
        
        // face-centered total energy density flux
        ctx.rhoEFlux[faceIdx] = 0.5 * ctx.faceArea[faceIdx] *
                                ((rhoELeft + ctx.pLeft[faceIdx]) * uBarLeft - uDotBLeft * bBarLeft +
                                (rhoERight + ctx.pRight[faceIdx]) * uBarRight - uDotBRight * bBarRight);
        
        // face-centered x-magnetic flux
        ctx.bXFlux[faceIdx] = 0.5 * ctx.faceArea[faceIdx] *
                              (ctx.uLeft[faceIdx] * bBarLeft - ctx.bXLeft[faceIdx] * uBarLeft +
                               ctx.uRight[faceIdx] * bBarRight - ctx.bXRight[faceIdx] * uBarRight);

        // face-centered y-magnetic flux
        ctx.bYFlux[faceIdx] = 0.5 * ctx.faceArea[faceIdx] *
                              (ctx.vLeft[faceIdx] * bBarLeft - ctx.bYLeft[faceIdx] * uBarLeft +
                               ctx.vRight[faceIdx] * bBarRight - ctx.bYRight[faceIdx] * uBarRight);

        // face-centered z-magnetic flux
        ctx.bZFlux[faceIdx] = 0.5 * ctx.faceArea[faceIdx] *
                              (ctx.wLeft[faceIdx] * bBarLeft - ctx.bZLeft[faceIdx] * uBarLeft +
                               ctx.wRight[faceIdx] * bBarRight - ctx.bZRight[faceIdx] * uBarRight);
    }

    FluxContext& ctx;
};

struct LowOrderGodunovKernel {
    LowOrderGodunovKernel(FluxContext& fluxContext) : ctx(fluxContext) {}

    inline void operator()(std::size_t const faceIdx) {
        // normal component of the velocity inside the cell
        auto uBar = ctx.uLeft[faceIdx] * ctx.faceNormalX[faceIdx] +
                    ctx.vLeft[faceIdx] * ctx.faceNormalY[faceIdx] +
                    ctx.wLeft[faceIdx] * ctx.faceNormalZ[faceIdx];

        // normal component of the magnetic field inside the cell
        auto bBar = ctx.bXLeft[faceIdx] * ctx.faceNormalX[faceIdx] +
                    ctx.bYLeft[faceIdx] * ctx.faceNormalY[faceIdx] +
                    ctx.bZLeft[faceIdx] * ctx.faceNormalZ[faceIdx];
        
        // cell-centered momentum densities inside the cell
        auto rhoU = ctx.rhoLeft[faceIdx] * ctx.uLeft[faceIdx];
        auto rhoV = ctx.rhoLeft[faceIdx] * ctx.vLeft[faceIdx];
        auto rhoW = ctx.rhoLeft[faceIdx] * ctx.wLeft[faceIdx];

        // cell-centered total energy density
        auto rhoE = ctx.rhoLeft[faceIdx] * (ctx.eLeft[faceIdx] + 0.5 * ctx.uuLeft[faceIdx]) + 0.5 * VACUUM_PERMITTIVITY_INV * ctx.bbLeft[faceIdx];

        // cell-centered dot(u,B)
        auto uDotB = ctx.uLeft[faceIdx] * ctx.bXLeft[faceIdx] +
                     ctx.vLeft[faceIdx] * ctx.bYLeft[faceIdx] +
                     ctx.wLeft[faceIdx] * ctx.bZLeft[faceIdx];

        // face-centered mass density flux
        ctx.rhoFlux[faceIdx] = ctx.faceArea[faceIdx] * ctx.rhoLeft[faceIdx] * uBar;

        // face-centered x-momentum density flux
        ctx.rhoUFlux[faceIdx] = ctx.faceArea[faceIdx] * (rhoU * uBar + ctx.pLeft[faceIdx] * ctx.faceNormalX[faceIdx] - ctx.bXLeft[faceIdx] * bBar);
                                                             
        // face-centered y-momentum density flux
        ctx.rhoVFlux[faceIdx] = ctx.faceArea[faceIdx] * (rhoV * uBar + ctx.pLeft[faceIdx] * ctx.faceNormalY[faceIdx] - ctx.bYLeft[faceIdx] * bBar);

        // face-centered z-momentum density flux
        ctx.rhoWFlux[faceIdx] = ctx.faceArea[faceIdx] * (rhoW * uBar + ctx.pLeft[faceIdx] * ctx.faceNormalZ[faceIdx] - ctx.bZLeft[faceIdx] * bBar);

        // face-centered total energy density flux
        ctx.rhoEFlux[faceIdx] = ctx.faceArea[faceIdx] * ((rhoE + ctx.pLeft[faceIdx]) * uBar - uDotB * bBar);
        
        // face-centered x-magnetic flux
        ctx.bXFlux[faceIdx] = ctx.faceArea[faceIdx] * (ctx.uLeft[faceIdx] * bBar - ctx.bXLeft[faceIdx] * uBar);

        // face-centered y-magnetic flux
        ctx.bYFlux[faceIdx] = ctx.faceArea[faceIdx] * (ctx.vLeft[faceIdx] * bBar - ctx.bYLeft[faceIdx] * uBar);

        // face-centered z-magnetic flux
        ctx.bZFlux[faceIdx] = ctx.faceArea[faceIdx] * (ctx.wLeft[faceIdx] * bBar - ctx.bZLeft[faceIdx] * uBar);
    }

    FluxContext& ctx;
};

class HLLC : public IFluxScheme {
public:
    void computeInterfaceFluxes(FluxContext& fluxContext) const {
        HLLCKernel kern(fluxContext);
    }
};

class KT : public IFluxScheme {
public:
    void computeInterfaceFluxes(FluxContext& fluxContext) const {
        KTKernel kern(fluxContext);
    }
};

class HighOrderGodunov : public IFluxScheme {
public:
    HighOrderGodunov(ExecutionController const& execCtrl) : execCtrl(execCtrl) {}
    void computeInterfaceFluxes(FluxContext& fluxContext) const {
        HighOrderGodunovKernel kern(fluxContext);
        
    }
    ExecutionController const& execCtrl;
};

class LowOrderGodunov : public IFluxScheme {
public:
    LowOrderGodunov(ExecutionController const& execCtrl) : execCtrl(execCtrl) {}
    void computeInterfaceFluxes(FluxContext& fluxContext) const {
        LowOrderGodunovKernel kern(fluxContext);
        // execCtrl.LaunchKernel(kern, fluxContext.faceArea.size());
    }
    ExecutionController const& execCtrl;
};

std::unique_ptr<IFluxScheme> flux_scheme_factory(Profile const& profile, ExecutionController const& execCtrl) {
    auto fluxSchemeName = profile.GetFluxScheme();
    switch (fluxSchemeName) {
    case FluxScheme::HLLC:
        return std::make_unique<HLLC>();
    case FluxScheme::KT:
        return std::make_unique<KT>();
    case FluxScheme::HOG:
        return std::make_unique<HighOrderGodunov>(execCtrl);
    case FluxScheme::LOG:
        return std::make_unique<LowOrderGodunov>(execCtrl);
    default:
        throw Error::INVALID_FLUX_SCHEME;
    }
}

} // namespace MHD
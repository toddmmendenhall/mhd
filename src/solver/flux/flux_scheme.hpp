#pragma once

#include <grid.hpp>
#include <memory>

namespace MHD {

class ExecutionController;
class Profile;
class ReconstructionContext;

struct FluxContext {
    FluxContext(IGrid const& grid, ReconstructionContext const& rc);

    std::size_t const numFaces;
    std::map<std::size_t, std::vector<std::size_t>> const& faceIdxToNodeIdxs;

    // properties of the faces
    std::vector<double> const& faceArea;
    std::vector<double> const& faceNormalX;
    std::vector<double> const& faceNormalY;
    std::vector<double> const& faceNormalZ;

    // Face-centered left states
    std::vector<double> const& rhoLeft;
    std::vector<double> const& uLeft;
    std::vector<double> const& vLeft;
    std::vector<double> const& wLeft;
    std::vector<double> const& pLeft;
    std::vector<double> const& eLeft;
    std::vector<double> const& csLeft;

    // Face-centered right states
    std::vector<double> const& rhoRight;
    std::vector<double> const& uRight;
    std::vector<double> const& vRight;
    std::vector<double> const& wRight;
    std::vector<double> const& pRight;
    std::vector<double> const& eRight;
    std::vector<double> const& csRight;

    // Face-centered fluxes
    std::vector<double> rhoFlux;
    std::vector<double> rhoUFlux;
    std::vector<double> rhoVFlux;
    std::vector<double> rhoWFlux;
    std::vector<double> rhoEFlux;
};

class IFlux {
public:
    virtual ~IFlux() = default;
    virtual void ComputeInterfaceFluxes(ExecutionController const& execCtrl) const = 0;
    FluxContext const& GetContext() const { return *m_context; }

protected:
    std::unique_ptr<FluxContext> m_context;
};

std::unique_ptr<IFlux> fluxFactory(Profile const& profile, IGrid const& grid, ReconstructionContext const& rc);

} // namespace MHD
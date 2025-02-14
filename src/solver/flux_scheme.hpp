#pragma once

#include <memory>

namespace MHD {

class Face;
class FluxContext;
class Profile;
class ExecutionController;

class IFluxScheme {
public:
    virtual ~IFluxScheme() = default;
    virtual void ComputeInterfaceFluxes(FluxContext& fluxContext) const = 0;
};

std::unique_ptr<IFluxScheme> fluxSchemeFactory(Profile const& profile, ExecutionController const& execCtrl);

} // namespace MHD
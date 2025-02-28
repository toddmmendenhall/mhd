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
    virtual void ComputeInterfaceFluxes(ExecutionController const& execCtrl) const = 0;
};

std::unique_ptr<IFluxScheme> fluxSchemeFactory(Profile const& profile, FluxContext& context);

} // namespace MHD
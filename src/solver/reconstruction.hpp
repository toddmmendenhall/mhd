#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
struct FluxContext;
class Profile;
struct ReconstructionContext;

class IReconstruction {
public:
    virtual ~IReconstruction() = default;
    virtual void compute(ReconstructionContext& ctx) const = 0;
};

std::unique_ptr<IReconstruction> reconstructionFactory(Profile const& profile, ExecutionController const& execCtrl);

} // namespace MHD
#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
class Profile;
struct ReconstructionContext;

class IReconstruction {
public:
    virtual ~IReconstruction() = default;
    virtual void Compute(ExecutionController const& execCtrl) = 0;
};

std::unique_ptr<IReconstruction> reconstructionFactory(Profile const& profile, ReconstructionContext& context);

} // namespace MHD
#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
class Profile;
struct ReconstructionContext;

class IReconstruction {
public:
    virtual ~IReconstruction() = default;
    virtual void ComputeReconstructedVariables(ExecutionController const& execCtrl, ReconstructionContext& context) const = 0;
};

std::unique_ptr<IReconstruction> reconstructionFactory(Profile const& profile);

} // namespace MHD
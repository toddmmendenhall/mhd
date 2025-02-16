#pragma once

#include <profile.hpp>

#include <memory>

namespace MHD {

class ExecutionController;
class BoundaryConditionContext;

class IBoundaryCondition {
public:
    virtual ~IBoundaryCondition() = default;
    virtual void Compute(ExecutionController const& execCtrl, BoundaryConditionContext& context) = 0;
};

std::unique_ptr<IBoundaryCondition> boundaryConditionFactory(Profile const& profile);

} // namespace MHD
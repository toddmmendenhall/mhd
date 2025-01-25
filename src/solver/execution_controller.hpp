#pragma once

#include <kernels.hpp>

namespace MHD
{

class ExecutionController
{
private:
    /* data */
public:
    ExecutionController(/* args */);
    ~ExecutionController();

    template <typename Func> void LaunchKernel(Func& kernel, size_t const n) const;
};

ExecutionController::ExecutionController(/* args */)
{
}

ExecutionController::~ExecutionController()
{
}

template <typename Func> void ExecutionController::LaunchKernel(Func& kernel, size_t const n) const
{
    for (size_t i = 0; i < n; ++i)
    {
        kernel(i);
    }
}

} // namespace MHD

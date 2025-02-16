#pragma once

#include <kernels.hpp>

#include <cstddef>

namespace MHD {

class ExecutionController {
public:
    ExecutionController();
    ~ExecutionController();

    template <typename Kernel> void LaunchKernel(Kernel& kernel, std::size_t const n) const;
    template <typename Kernel> void LaunchKernel(Kernel& kernel, std::size_t const m, std::size_t const n) const;
};

} // namespace MHD

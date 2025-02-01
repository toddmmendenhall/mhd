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

template <typename Kernel> void ExecutionController::LaunchKernel(Kernel& kernel, std::size_t const n) const {
    for (std::size_t i = 0; i < n; ++i) {
        kernel(i);
    }
}

template <typename Kernel> void ExecutionController::LaunchKernel(Kernel& kernel, std::size_t const m, std::size_t const n) const {
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            kernel(i, j);
        }
    }
}

} // namespace MHD

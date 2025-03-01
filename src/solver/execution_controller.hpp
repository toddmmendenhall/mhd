#pragma once

#include <cstddef>
#include <vector>

namespace MHD {

class ExecutionController {
public:
    template <typename Kernel> void LaunchKernel(Kernel& kernel, std::size_t const n) const {
        for (std::size_t i = 0; i < n; ++i) {
            kernel(i);
        }
    }
    
    template <typename Kernel> void LaunchKernel(Kernel& kernel, std::size_t const m, std::size_t const n) const {
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                kernel(i, j);
            }
        }
    }

    template <typename Kernel> void LaunchKernel(Kernel& kernel, std::vector<std::size_t> const& idxs) const {
        for (std::size_t i : idxs) {
            kernel(i);
        }
    }
};

} // namespace MHD

#pragma once

#include <vector>

namespace MHD {

struct CentralFiniteDifference1DKernel {
    CentralFiniteDifference1DKernel(std::vector<double> const& rho, double const dxInverse, std::vector<double>& drhodx) :
        rho(rho), dxInverse(dxInverse), drhodx(drhodx) {}
    
    void operator()(std::size_t idx) {
        drhodx[idx] = (rho[idx+1] - rho[idx-1]) * 0.5 * dxInverse;
    }

    std::vector<double> const& rho;
    double const dxInverse;
    std::vector<double>& drhodx;
};

} // namespace MHD
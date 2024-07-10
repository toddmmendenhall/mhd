#include "cartesian.hpp"

#include "profile.hpp"

#include <memory>
#include <vector>

namespace MHD {

CartesianGrid::CartesianGrid(std::unique_ptr<GridProfile> const& gridProfile) {
    auto d = gridProfile->m_dimension;
    auto l = gridProfile->m_limits;
    auto s = gridProfile->m_spacing;
    std::vector<std::vector<double>> values;
    for (std::size_t i = 0; i < d; ++i) {
        std::vector<double> dimensionValues;
        for (double value = l[2 * i]; value <= l[2 * i + 1]; value += s[i]) {
            dimensionValues.push_back(value);
        }
        values.push_back(dimensionValues);
    }
}

CartesianGrid::~CartesianGrid() {}

} // namespace MHD
#include <1d.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>
#include <vector>

namespace MHD {

Cartesian1DGrid::Cartesian1DGrid(Profile const& profile) : IGrid(profile) {
    auto& bounds = profile.GetGridBounds();
    auto& spacing = profile.GetGridSpacing();
    auto numXNodes = (bounds[1] - bounds[0]) / spacing[0];
    
    // Append interior nodes
    for (std::size_t i = 0; i < numXNodes; ++i) {
        m_nodePositions.push_back(bounds[0] + (i + 0.5) * spacing[0]);
    }

    // Insert boundary nodes
    m_nodePositions.insert(m_nodePositions.begin(), bounds[0] - 0.5 * spacing[0]);
    m_nodePositions.insert(m_nodePositions.end(), bounds[1] + 0.5 * spacing[0]);
}

Cartesian1DGrid::~Cartesian1DGrid() = default;

void Cartesian1DGrid::SomeMethod() {}

} // namespace MHD

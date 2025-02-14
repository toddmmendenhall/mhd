#include <1d.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>
#include <vector>

namespace MHD {

Cartesian1DGrid::Cartesian1DGrid(Profile const& profile) {
    auto& bounds = profile.m_gridBoundsOption;
    auto& spacing = profile.m_gridSpacingsOption;

    auto numInteriorNodes = (bounds[1] - bounds[0]) / spacing[0];

    // Append interior nodes
    for (std::size_t i = 0; i < numInteriorNodes; ++i) {
        m_nodes.push_back({bounds[0] + (i + 0.5) * spacing[0], 0.0, 0.0});
    }

    // Insert boundary nodes
    m_nodes.insert(m_nodes.begin(), {bounds[0] - 0.5 * spacing[0], 0.0, 0.0});
    m_nodes.insert(m_nodes.end(), {bounds[1] + 0.5 * spacing[0], 0.0, 0.0});

    // Each face has a "left" and "right" node, so we store those indices
    for (std::size_t i = 0; i < m_nodes.size() - 1; ++i) {
        m_faceToNodeIndices.push_back({i, i+1});
    }

    // The area of each face is determined from the other two dimensions
    double const faceArea = spacing[1] * spacing[2];
    for (std::size_t i = 0; i < m_faceToNodeIndices.size(); ++i) {
        m_faceAreas.push_back(faceArea);
    }

    // The domain is assumed to lie on the x-axis
    for (std::size_t i = 0; i < m_faceToNodeIndices.size(); ++i) {
        m_faceNormalsX.push_back(1.0);
        m_faceNormalsY.push_back(0.0);
        m_faceNormalsZ.push_back(0.0);
    }
}

} // namespace MHD

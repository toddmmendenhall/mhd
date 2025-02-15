#include <1d.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>
#include <vector>

namespace MHD {

Cartesian1DGrid::Cartesian1DGrid(Profile const& profile) {
    auto& bounds = profile.m_gridBoundsOption;
    auto& spacing = profile.m_gridSpacingsOption;

    std::size_t numInteriorNodes = (bounds[1] - bounds[0]) / spacing[0];

    // Append cell-centered interior nodes
    for (std::size_t i = 0; i < numInteriorNodes; ++i) {
        m_nodes.push_back({bounds[0] + (i + 0.5) * spacing[0], 0.0, 0.0});
    }

    // Append face-centered boundary nodes
    m_nodes.push_back({bounds[0], 0.0, 0.0});
    m_nodes.push_back({bounds[1], 0.0, 0.0});

    // Each face has a "left" and "right" node, so we store those indices
    for (std::size_t i = 0; i < numInteriorNodes + 1; ++i) {
        if (i == 0) {
            // The first interior node has the lower boundary on the left
            m_faceToNodeIndices.push_back({numInteriorNodes, i});
        } else if (i == numInteriorNodes - 1) {
            // The last interior node has the upper boundary on the right
            m_faceToNodeIndices.push_back({i - 1, numInteriorNodes + 1});
        } else {
            m_faceToNodeIndices.push_back({i - 1, i});
        }
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

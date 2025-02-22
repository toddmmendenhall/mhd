#include <1d.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>
#include <vector>

namespace MHD {

Cartesian1DGrid::Cartesian1DGrid(Profile const& profile) {
    auto& bounds = profile.m_gridBoundsOption;
    cellSize = profile.m_gridSpacingsOption;

    numInteriorNodes = (bounds[1] - bounds[0]) / cellSize[0];

    // Append cell-centered interior nodes
    for (std::size_t i = 0; i < numInteriorNodes; ++i) {
        m_nodes.push_back({bounds[0] + (i + 0.5) * cellSize[0], 0.0, 0.0});
    }

    // Append symmetric-across-the-boundary fictitious boundary nodes
    m_nodes.push_back({bounds[0] - 0.5 * cellSize[0], 0.0, 0.0});
    m_nodes.push_back({bounds[1] + 0.5 * cellSize[0], 0.0, 0.0});

    // Each face has a "left" and "right" node, so we store those indices
    for (std::size_t i = 0; i < m_nodes.size() - 1; ++i) {
        if (i == 0) {
            // The left boundary face
            m_faceToNodeIndices.push_back({numInteriorNodes, i});
        } else if (i == numInteriorNodes) {
            // The right boundary face
            m_faceToNodeIndices.push_back({i - 1, numInteriorNodes + 1});
        } else {
            // Interior faces
            m_faceToNodeIndices.push_back({i - 1, i});
        }
    }

    // The area of each face is determined from the other two dimensions
    double const faceArea = cellSize[1] * cellSize[2];
    for (std::size_t i = 0; i < m_faceToNodeIndices.size(); ++i) {
        m_faceAreas.push_back(faceArea);
    }

    // The domain is assumed to lie on the x-axis
    for (std::size_t i = 0; i < m_faceToNodeIndices.size(); ++i) {
        if (i == 0) {
            // The left boundary normal faces out of the domain
            m_faceNormalsX.push_back(-1.0);
        } else if (i == numInteriorNodes) {
            // The right boundary normal faces out of the domain
            m_faceNormalsX.push_back(1.0);
        } else {
            m_faceAreas.push_back(1.0);
        }
        m_faceNormalsY.push_back(0.0);
        m_faceNormalsZ.push_back(0.0);
    }

    for (std::size_t i = 0; i < numInteriorNodes; ++i) {
        cellToFaceIndices.push_back({i, i + 1});
    }
}

} // namespace MHD

#include <1d.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>
#include <vector>

namespace MHD {

Cartesian1DGrid::Cartesian1DGrid(Profile const& profile) {
    auto& bounds = profile.m_gridBoundsOption;
    cellSize = profile.m_gridSpacingsOption;
    numCells = (bounds[1] - bounds[0]) / cellSize[0];
    numFaces = numCells + 1;

    // Append cell-centered real nodes
    for (std::size_t i = 0; i < numCells; ++i) {
        m_nodes.push_back({bounds[0] + (i + 0.5) * cellSize[0], 0.0, 0.0});
    }

    // Add cell-centered ghost nodes at the ends of the domain
    m_nodes.insert(m_nodes.begin(), {bounds[0] - 0.5 * cellSize[0], 0.0, 0.0});
    m_nodes.push_back({bounds[1] + 0.5 * cellSize[0], 0.0, 0.0});

    startIdx = 1;

    // Each face has a "left" and "right" node, so we store those indices
    for (FaceIdx i = 0; i < numFaces; ++i) {
        NodeIdxs nodeIdxs;
        nodeIdxs.left = i;
        nodeIdxs.right = i + 1;
        nodeIdxs.leftMinusOne = nodeIdxs.left - 1;
        nodeIdxs.rightPlusOne = nodeIdxs.right + 1;
        nodeIdxs.inner = -1;
        nodeIdxs.outer = -1;
        if (i == 0) {
            nodeIdxs.isBoundary = true;
            nodeIdxs.inner = i + 1;
            nodeIdxs.outer = i;
            nodeIdxs.leftMinusOne = i;
        } else if (i == numFaces - 1) {
            nodeIdxs.isBoundary = true;
            nodeIdxs.inner = i;
            nodeIdxs.outer = i + 1;
            nodeIdxs.rightPlusOne = i;
        }
        m_faceIdxToNodeIdxs[i] = nodeIdxs;
    }

    // The area of each face is determined from the other two dimensions
    double const faceArea = cellSize[1] * cellSize[2];
    for (std::size_t i = 0; i < numFaces; ++i) {
        m_faceAreas.push_back(faceArea);
    }

    // The domain is assumed to lie on the x-axis
    for (std::size_t i = 0; i < numFaces; ++i) {
        if (i == 0) {
            // The left boundary normal faces out of the domain
            m_faceNormalsX.push_back(-1.0);
        } else if (i == numFaces - 1) {
            // The right boundary normal faces out of the domain
            m_faceNormalsX.push_back(1.0);
        } else {
            m_faceNormalsX.push_back(1.0);
        }
        m_faceNormalsY.push_back(0.0);
        m_faceNormalsZ.push_back(0.0);
    }

    for (std::size_t i = 0; i < numCells; ++i) {
        cellToFaceIndices.push_back({i, i + 1});
    }
}

} // namespace MHD

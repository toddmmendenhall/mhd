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

    // Append cell-centered nodes for each cell
    for (std::size_t i = 0; i < numCells; ++i) {
        m_nodes.push_back({bounds[0] + (i + 0.5) * cellSize[0], 0.0, 0.0});
    }

    // Insert nodes for left ghost cells
    // m_nodes.insert(m_nodes.begin(), {bounds[0] - 0.5 * cellSize[0], 0.0, 0.0}); // 0
    // m_nodes.insert(m_nodes.begin(), {bounds[0] - 1.5 * cellSize[0], 0.0, 0.0}); // 1
    // m_nodes.insert(m_nodes.begin(), {bounds[0] - 2.5 * cellSize[0], 0.0, 0.0}); // 2
    m_nodes.push_back({bounds[0] - 0.5 * cellSize[0], 0.0, 0.0}); // 10
    m_nodes.push_back({bounds[0] - 1.5 * cellSize[0], 0.0, 0.0}); // 11
    m_nodes.push_back({bounds[0] - 2.5 * cellSize[0], 0.0, 0.0}); // 12

    // Append nodes for right ghost cells
    m_nodes.push_back({bounds[1] + 0.5 * cellSize[0], 0.0, 0.0}); // 13
    m_nodes.push_back({bounds[1] + 1.5 * cellSize[0], 0.0, 0.0}); // 14
    m_nodes.push_back({bounds[1] + 2.5 * cellSize[0], 0.0, 0.0}); // 15

    // Each face has a "left" and "right" node, so we store those indices
    for (std::size_t i = 0; i < numFaces; ++i) {
        if (i == 0) {
            // The left boundary face
            m_faceToNodeIndices.push_back({numCells+1, numCells, 0, 1});
            boundaryFaceToBoundaryCellIndices.push_back(numCells);
            boundaryFaceToInteriorCellIndices.push_back(i);
        } else if (i == 1) {
            m_faceToNodeIndices.push_back({numCells, 0, 1, 2});
        } else if (i == numFaces - 2) {
            m_faceToNodeIndices.push_back({i-2, i-1, i, i+4});
        } else if (i == numFaces - 1) {
            // The right boundary face
            m_faceToNodeIndices.push_back({i-2, i-1, i+3, i+4});
            boundaryFaceToBoundaryCellIndices.push_back(numCells + 3);
            boundaryFaceToInteriorCellIndices.push_back(i - 1);
        } else {
            // Interior faces
            m_faceToNodeIndices.push_back({i-2, i-1, i, i+1});
        }
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

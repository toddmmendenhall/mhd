#include <1d.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>
#include <vector>

namespace MHD {

Cartesian1DGrid::Cartesian1DGrid(Profile const& profile) {
    auto& bounds = profile.m_gridBoundsOption;
    m_cellSize = profile.m_gridSpacingsOption;
    m_numCells = (bounds[1] - bounds[0]) / m_cellSize[0];
    m_numFaces = m_numCells + 1;
    m_numBoundaries = 2;

    // Internal nodes correspond to the cell centers
    for (std::size_t i = 0; i < m_numCells; ++i) {
        m_nodes.push_back({bounds[0] + (i + 0.5) * m_cellSize[0], 0.0, 0.0});
    }
    
    // External nodes correspond to the ghost cells
    m_nodes.push_back({bounds[0] - 0.5 * m_cellSize[0], 0.0, 0.0});
    m_nodes.push_back({bounds[1] + 0.5 * m_cellSize[0], 0.0, 0.0});

    // Each face has a "left" and "right" cell
    // Boundary faces have an "inner" and "outer" cell
    for (std::size_t i = 0; i < m_numFaces; ++i) {
        m_faceIdxs.push_back(i);
        if (i == 0) {
            m_faceIdxToCellIdxs[i] = {m_numCells, i, m_numCells, i + 1};
            m_boundaryIdxs.push_back(i);
            m_boundaryIdxToCellIdxs[i] = {i, m_numCells};
        } else if (i == 1) {
            m_faceIdxToCellIdxs[i] = {i - 1, i, m_numCells, i + 1};
        } else if (i == m_numFaces - 2) {
            m_faceIdxToCellIdxs[i] = {i - 1, i, i - 2, m_numCells + 1};
        } else if (i == m_numFaces - 1) {
            m_faceIdxToCellIdxs[i] = {i - 1, m_numCells + 1, i - 2, m_numCells + 1};
            m_boundaryIdxs.push_back(i);
            m_boundaryIdxToCellIdxs[i] = {i - 1, m_numCells + 1};
        } else {
            m_faceIdxToCellIdxs[i] = {i - 1, i, i - 2, i + 1};
        }
    }

    //  Each cell has a "left" and "right" face
    for (std::size_t i = 0; i < m_numCells; ++i) {
        m_cellIdxToFaceIdxs[i] = {i, i + 1};
    }

    // The area of each face is determined from the other two dimensions
    double const faceArea = 1.0;

    // The domain is assumed to lie on the x-axis
    for (std::size_t i = 0; i < m_numFaces; ++i) {
        m_faceAreas.push_back(faceArea);
        m_faceNormalsX.push_back(1.0);
        m_faceNormalsY.push_back(0.0);
        m_faceNormalsZ.push_back(0.0);
    }
}

} // namespace MHD

#pragma once

#include <map>
#include <memory>
#include <vector>

namespace MHD {

class Profile;

class IGrid {
public:
    virtual ~IGrid() = default;

    std::vector<std::array<double, 3>> const& Nodes() const { return m_nodes; }
    std::size_t const NumCells() const { return m_numCells; }
    std::size_t const NumFaces() const { return m_numFaces; }
    std::size_t const NumBoundaries() const { return m_numBoundaries; }
    std::size_t const NumNodes() const { return m_numCells + m_numBoundaries; }
    std::map<std::size_t, std::vector<std::size_t>> const& FaceIdxToCellIdxs() const { return m_faceIdxToCellIdxs; }
    std::map<std::size_t, std::vector<std::size_t>> const& CellIdxToFaceIdxs() const { return m_cellIdxToFaceIdxs; }
    std::map<std::size_t, std::vector<std::size_t>> const& BoundaryIdxToCellIdxs() const { return m_boundaryIdxToCellIdxs; }
    std::vector<std::size_t> const& BoundaryIdxs() const { return m_boundaryIdxs; }
    std::vector<std::size_t> const& FaceIdxs() const { return m_faceIdxs; }
    std::vector<double> const& FaceAreas() const { return m_faceAreas; }
    std::vector<double> const& FaceNormalX() const { return m_faceNormalsX; }
    std::vector<double> const& FaceNormalY() const { return m_faceNormalsY; }
    std::vector<double> const& FaceNormalZ() const { return m_faceNormalsZ; }

    std::vector<double> const& CellSize() const { return m_cellSize; }

protected:
    std::vector<std::array<double, 3>> m_nodes;
    std::size_t m_numCells;
    std::size_t m_numFaces;
    std::size_t m_numBoundaries;
    std::map<std::size_t, std::vector<std::size_t>> m_faceIdxToCellIdxs;
    std::map<std::size_t, std::vector<std::size_t>> m_cellIdxToFaceIdxs;
    std::map<std::size_t, std::vector<std::size_t>> m_boundaryIdxToCellIdxs;
    std::vector<std::size_t> m_boundaryIdxs;
    std::vector<std::size_t> m_faceIdxs;
    std::vector<double> m_faceAreas;
    std::vector<double> m_faceNormalsX;
    std::vector<double> m_faceNormalsY;
    std::vector<double> m_faceNormalsZ;

    std::vector<double> m_cellSize;
};

std::unique_ptr<IGrid> gridFactory(Profile const& profile);

} // namespace MHD

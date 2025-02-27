#pragma once

#include <memory>
#include <vector>

namespace MHD {

class Profile;

class IGrid {
public:
    virtual ~IGrid() = default;

    std::vector<std::array<double, 3>> const& Nodes() const { return m_nodes; }
    std::size_t const NumCells() const { return numCells; }
    std::size_t const NumNodes() const { return m_nodes.size(); }
    std::vector<std::array<std::size_t, 2>> const& FaceToNodeIndices() const { return m_faceToNodeIndices; }
    std::vector<std::array<std::size_t, 2>> const& CellToFaceIndices() const { return cellToFaceIndices; }
    std::vector<std::size_t> const& BoundaryFaceToBoundaryCellIndices() const {
        return boundaryFaceToBoundaryCellIndices;
    }
    std::vector<std::size_t> const& BoundaryFaceToInteriorCellIndices() const {
        return boundaryFaceToInteriorCellIndices;
    }
    std::vector<double> const& FaceAreas() const { return m_faceAreas; }
    std::vector<double> const& FaceNormalX() const { return m_faceNormalsX; }
    std::vector<double> const& FaceNormalY() const { return m_faceNormalsY; }
    std::vector<double> const& FaceNormalZ() const { return m_faceNormalsZ; }
    std::size_t const NumFaces() const { return numFaces; }
    std::vector<double> const& CellSize() const { return cellSize; }

protected:
    std::vector<std::array<double, 3>> m_nodes;
    std::size_t numCells;
    std::size_t numFaces;
    std::vector<std::array<std::size_t, 2>> m_faceToNodeIndices;
    std::vector<std::array<std::size_t, 2>> cellToFaceIndices;
    std::vector<std::size_t> boundaryFaceToBoundaryCellIndices;
    std::vector<std::size_t> boundaryFaceToInteriorCellIndices;
    std::vector<double> m_faceAreas;
    std::vector<double> m_faceNormalsX;
    std::vector<double> m_faceNormalsY;
    std::vector<double> m_faceNormalsZ;
    std::vector<double> cellSize;
};

std::unique_ptr<IGrid> gridFactory(Profile const& profile);

} // namespace MHD

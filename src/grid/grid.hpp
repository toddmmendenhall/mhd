#pragma once

#include <map>
#include <memory>
#include <vector>

namespace MHD {

class Profile;


using FaceIdx = std::size_t;
using NodeIdx = std::size_t;
using CellIdx = std::size_t;
struct NodeIdxs {
    NodeIdx left;
    NodeIdx right;
    NodeIdx leftMinusOne;
    NodeIdx rightPlusOne;
    bool isBoundary = false;
    NodeIdx inner;
    NodeIdx outer;
};
using FaceIdxToNodeIdxs = std::map<FaceIdx, NodeIdxs>;
struct FaceIdxs {
    FaceIdx left;
    FaceIdx right;
};

class IGrid {
public:
    virtual ~IGrid() = default;

    std::vector<std::array<double, 3>> const& Nodes() const { return m_nodes; }
    std::size_t const NumCells() const { return numCells; }
    std::size_t const NumNodes() const { return m_nodes.size(); }
    FaceIdxToNodeIdxs const& GetFaceIdxToNodeIdxs() const { return m_faceIdxToNodeIdxs; }
    std::map<CellIdx, FaceIdxs> const& GetCellIdxToFaceIdxs() const { return cellIdxToFaceIdxs; }
    std::vector<std::size_t> const& GetCellIdxs() const { return cellIdxs; }

    std::vector<double> const& FaceAreas() const { return m_faceAreas; }
    std::vector<double> const& FaceNormalX() const { return m_faceNormalsX; }
    std::vector<double> const& FaceNormalY() const { return m_faceNormalsY; }
    std::vector<double> const& FaceNormalZ() const { return m_faceNormalsZ; }
    std::size_t const NumFaces() const { return numFaces; }
    std::vector<double> const& CellSize() const { return cellSize; }

protected:
    std::vector<std::size_t> cellIdxs;
    std::vector<std::array<double, 3>> m_nodes;
    std::size_t numCells;
    std::size_t numFaces;
    FaceIdxToNodeIdxs m_faceIdxToNodeIdxs;
    std::map<CellIdx, FaceIdxs> cellIdxToFaceIdxs;
    std::vector<double> m_faceAreas;
    std::vector<double> m_faceNormalsX;
    std::vector<double> m_faceNormalsY;
    std::vector<double> m_faceNormalsZ;
    std::vector<double> cellSize;
};

std::unique_ptr<IGrid> gridFactory(Profile const& profile);

} // namespace MHD

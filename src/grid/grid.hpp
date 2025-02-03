#pragma once

#include <memory>
#include <vector>

namespace MHD {

class Profile;

class IGrid {
public:
    virtual ~IGrid() = default;

    std::vector<std::array<double, 3>> const& Nodes() const { return m_nodes; }
    std::vector<std::array<std::size_t, 2>> const& FaceToNodeIndices() const { return m_faceToNodeIndices; }
    std::vector<double> const& FaceAreas() const { return m_faceAreas; }
    std::vector<double> const& FaceNormalX() const { return m_faceNormalsX; }
    std::vector<double> const& FaceNormalY() const { return m_faceNormalsY; }
    std::vector<double> const& FaceNormalZ() const { return m_faceNormalsZ; }
    std::size_t const NumFaces() const { return m_faceToNodeIndices.size(); }

protected:
    std::vector<std::array<double, 3>> m_nodes;
    std::vector<std::array<std::size_t, 2>> m_faceToNodeIndices;
    std::vector<double> m_faceAreas;
    std::vector<double> m_faceNormalsX;
    std::vector<double> m_faceNormalsY;
    std::vector<double> m_faceNormalsZ;
};

std::unique_ptr<IGrid> GridFactory(Profile const& profile);

} // namespace MHD

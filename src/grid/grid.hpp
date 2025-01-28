#pragma once

#include <memory>
#include <vector>

namespace MHD {

class Profile;

class IGrid {
public:
    IGrid(Profile const& profile) { m_numFaces = 10; }
    ~IGrid() = default;

    virtual double const GetDx() const = 0;
    virtual std::vector<double> const& GetNodePositions() const = 0;
    std::size_t const NumFaces() const { return m_numFaces; }

    virtual void SomeMethod() = 0;

protected:
    std::size_t m_numFaces;
};

std::unique_ptr<IGrid> grid_factory(Profile const& profile);

} // namespace MHD

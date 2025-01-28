#pragma once

#include <grid.hpp>

#include <memory>
#include <vector>

namespace MHD {

class Profile;

class Cartesian1DGrid : public IGrid {
public:
    Cartesian1DGrid(Profile const& profile);
    ~Cartesian1DGrid();
    double const GetDx() const {return m_dx;}
    std::vector<double> const& GetNodePositions() const {return m_nodePositions;}

    void SomeMethod() override;

private:
    std::vector<double> m_nodePositions;
    double m_dx;
};

} // namespace MHD

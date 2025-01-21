#pragma once

#include <memory>
#include <vector>

namespace MHD {

class Profile;

class Grid {
public:
    Grid(Profile const& profile);
    ~Grid();

    virtual double const GetDx() const = 0;
    virtual std::vector<double> const& GetNodePositions() const = 0;

    virtual void SomeMethod() = 0;
};

} // namespace MHD

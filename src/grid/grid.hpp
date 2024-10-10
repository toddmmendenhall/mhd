#pragma once

#include <memory>

namespace MHD {

class Profile;

class Grid {
public:
    Grid(Profile const& profile);
    ~Grid();

    virtual void SomeMethod() = 0;
};

} // namespace MHD

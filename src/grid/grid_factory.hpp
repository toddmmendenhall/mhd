#pragma once

#include <memory>

namespace MHD {

class Grid;
class Profile;

class GridFactory {
public:
    GridFactory();
    ~GridFactory();

    std::unique_ptr<Grid> CreateGrid(Profile const& profile) const;
};

} // namespace MHD

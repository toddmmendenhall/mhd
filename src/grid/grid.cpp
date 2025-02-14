#include <cartesian/1d.hpp>
#include <error.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>

namespace MHD {

std::unique_ptr<IGrid> GridFactory(Profile const& profile) {
    auto geometry = profile.m_gridGeometryOption;
    auto dimension = profile.m_gridDimensionOption;

    if (geometry == Geometry::CARTESIAN) {
        if (dimension == Dimension::ONE) {
            return std::make_unique<Cartesian1DGrid>(profile);
        } else {
            throw Error::INVALID_GRID_DIMENSION;
        }
    } else {
        throw Error::INVALID_GRID_GEOMETRY;
    }
}

} // namespace MHD

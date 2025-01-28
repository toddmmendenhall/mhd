#include <cartesian/1d.hpp>
// #include <cartesian/2d.hpp>
// #include <cartesian/3d.hpp>
#include <error.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>

namespace MHD {

std::unique_ptr<IGrid> grid_factory(Profile const& profile) {
    Dimension dimension = profile.GetGridDimension();
    Geometry geometry = profile.GetGridGeometry();

    if (geometry == Geometry::CARTESIAN) {
        if (dimension == Dimension::ONE) {
            return std::make_unique<Cartesian1DGrid>(profile);
        // } else if (dimension == Dimension::TWO) {
            // return std::make_unique<Cartesian2DGrid>(profile);
        // } else if (dimension == Dimension::THREE) {
            // return std::make_unique<Cartesian3DGrid>(profile);
        } else {
            throw Error::INVALID_GRID_DIMENSION;
        }
    } else {
        throw Error::INVALID_GRID_GEOMETRY;
    }
}

} // namespace MHD

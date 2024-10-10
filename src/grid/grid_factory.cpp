#include <1d.hpp>
#include <2d.hpp>
#include <3d.hpp>
#include <grid.hpp>
#include <grid_factory.hpp>
#include <profile.hpp>
#include <profile_options.hpp>

#include <memory>
#include <vector>

namespace MHD {

GridFactory::GridFactory() = default;

GridFactory::~GridFactory() = default;

std::unique_ptr<Grid> GridFactory::CreateGrid(Profile const& profile) const {
    Dimension dimension = profile.GetGridDimension();
    Geometry geometry = profile.GetGridGeometry();

    if (geometry == Geometry::CARTESIAN) {
        if (dimension == Dimension::ONE) {
            return std::make_unique<Cartesian1DGrid>(profile);
        }
        if (dimension == Dimension::TWO) {
            return std::make_unique<Cartesian2DGrid>(profile);
        }
        if (dimension == Dimension::THREE) {
            return std::make_unique<Cartesian3DGrid>(profile);
        }
    }
    return nullptr;
}

} // namespace MHD

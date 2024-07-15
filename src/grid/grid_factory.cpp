#include "cartesian.hpp"
#include "grid.hpp"
#include "grid_factory.hpp"
#include "profile.hpp"
#include "profile_options.hpp"

#include <vector>

namespace MHD {

Grid* GridFactory::CreateGrid(Profile* profile) const {
    Dimension dimension = profile->GetGridDimension();
    Geometry geometry = profile->GetGridGeometry();

    if (dimension == Dimension::ONE && geometry == Geometry::CARTESIAN) {
        return new Cartesian1DGrid(profile);
    }
    if (dimension == Dimension::TWO && geometry == Geometry::CARTESIAN) {
        return new Cartesian2DGrid(profile);
    }
    if (dimension == Dimension::THREE && geometry == Geometry::CARTESIAN) {
        return new Cartesian3DGrid(profile);
    }

    return nullptr;
}

} // namespace MHD

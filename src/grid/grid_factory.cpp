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
        return new Cartesian1DGrid(profile->GetGridBounds(), profile->GetGridSpacing(), profile->GetGridBoundaryConditions());
    }
    if (dimension == Dimension::TWO && geometry == Geometry::CARTESIAN) {
        return new Cartesian2DGrid(profile->GetGridBounds(), profile->GetGridSpacing(), profile->GetGridBoundaryConditions());
    }
    if (dimension == Dimension::THREE && geometry == Geometry::CARTESIAN) {
        return new Cartesian3DGrid(profile->GetGridBounds(), profile->GetGridSpacing(), profile->GetGridBoundaryConditions());
    }

    return nullptr;
}

} // namespace MHD
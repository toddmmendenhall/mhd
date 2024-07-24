#include "1d.hpp"
#include "2d.hpp"
#include "3d.hpp"
#include "grid_factory.hpp"
#include "grid_impl.hpp"
#include "profile.hpp"
#include "profile_options.hpp"

#include <memory>
#include <vector>

namespace MHD {

std::unique_ptr<GridImpl> GridFactory::CreateGrid(std::unique_ptr<Profile> const& profile) const {
    Dimension dimension = profile->GetGridDimension();
    Geometry geometry = profile->GetGridGeometry();

    if (dimension == Dimension::ONE && geometry == Geometry::CARTESIAN) {
        return std::make_unique<Cartesian1DGrid>(profile);
    }

    if (dimension == Dimension::TWO && geometry == Geometry::CARTESIAN) {
        return std::make_unique<Cartesian2DGrid>(profile);
    }

    if (dimension == Dimension::THREE && geometry == Geometry::CARTESIAN) {
        return std::make_unique<Cartesian3DGrid>(profile);
    }

    return nullptr;
}

} // namespace MHD

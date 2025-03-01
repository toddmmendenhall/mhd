#include <1d.hpp>
#include <error.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>

namespace MHD {

std::unique_ptr<IGrid> gridFactory(Profile const& profile) {
    if (Dimension::ONE == profile.m_gridDimensionOption) {
        return std::make_unique<Cartesian1DGrid>(profile);
    }
    return nullptr;
}

} // namespace MHD

#include "domain_impl.hpp"
#include "profile.hpp"
#include "profile_options.hpp"

namespace MHD {

DomainImpl::DomainImpl(Profile const* profile) {
    profile->GetGridDimension();
    profile->GetGridGeometry();
    profile->GetGridBounds();
    profile->GetGridSpacing();
    profile->GetGridBoundaryConditions();
}

} // namespace MHD

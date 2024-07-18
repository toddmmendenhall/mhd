#include "domain.hpp"
#include "domain_impl.hpp"
#include "profile.hpp"

namespace MHD {

Domain::Domain(Profile const* profile) {
    m_domain = new DomainImpl(profile);
}

Domain::~Domain() {
    delete m_domain;
}

} // namespace MHD

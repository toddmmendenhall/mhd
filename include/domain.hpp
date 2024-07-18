#pragma once

#include "profile.hpp"

namespace MHD {

class DomainImpl;

class Domain {
public:
    Domain(Profile const* profile);
    ~Domain();

private:
   DomainImpl* m_domain;
};

} // namespace MHD
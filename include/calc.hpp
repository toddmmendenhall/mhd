#pragma once

#include "domain.hpp"
#include "grid.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

class Calc {
public:
    Calc();
    ~Calc();

    void Run();

    Profile* GetProfile() const;
    Grid* GetGrid() const;    

private:
    Profile* m_profile;
    Grid* m_grid;
    std::unique_ptr<Domain> m_domain;
};

} // namespace MHD

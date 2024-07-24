#include "calc.hpp"
#include "grid.hpp"
#include "profile.hpp"

#include <iostream>
#include <memory>

namespace MHD {

Calc::Calc() {
    m_profile = std::make_unique<Profile>();
}

void Calc::Run() {
    SetupCalc();
    std::cout << "Running calc..." << std::endl;
}

void Calc::SetupCalc() {
    m_grid = std::make_unique<Grid>(m_profile);
}

std::unique_ptr<Profile> const& Calc::GetProfile() const {
    return m_profile;
}

std::unique_ptr<Grid> const& Calc::GetGrid() const {
    return m_grid;
}

} // namespace MHD

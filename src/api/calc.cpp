#include <calc.hpp>
#include <grid.hpp>
#include <grid_factory.hpp>
#include <profile.hpp>

#include <memory>

namespace MHD {

Calc::Calc() {
    m_profile = std::make_unique<Profile>();

    GridFactory gridFactory;
    m_grid = gridFactory.CreateGrid(*m_profile);
}

Calc::~Calc() = default;

void Calc::Run() {}

} // namespace MHD

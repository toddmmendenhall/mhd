#include "api_safe.hpp"
#include "calc.hpp"
#include "grid.hpp"
#include "logger.hpp"
#include "profile.hpp"

#include <memory>

Calc::Calc() {
    m_logger = std::make_unique<Logger>();
    m_profile = std::make_unique<Profile>();
    m_grid = std::make_unique<Grid>(m_profile);
}

Calc::~Calc() {}

void Calc::Run() {
    api_safe(m_logger, [&]() {
        // do stuff here
    });
}

#pragma once

#include "grid.hpp"
#include "logger.hpp"
#include "profile.hpp"

#include <memory>

class Calc {
public:
    Calc();
    ~Calc();

    void Run();

private:
    std::unique_ptr<Logger> m_logger;
    std::unique_ptr<Profile> m_profile;
    std::unique_ptr<Grid> m_grid;
};

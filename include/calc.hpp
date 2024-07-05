#pragma once

#include "error.hpp"

#include <memory>

class Logger;

class Calc {
public:
    Calc();
    ~Calc();

    Error Run();

private:
    std::unique_ptr<Logger> m_logger;
};

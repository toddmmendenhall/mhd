#pragma once

#include "error.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

class Logger;

class Calc {
public:
    Calc();
    ~Calc();

    Error Run();

private:
    std::unique_ptr<Logger> m_logger;
    std::unique_ptr<Profile> m_profile;
};

}

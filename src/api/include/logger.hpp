#pragma once

#include "error.hpp"
#include "exception.hpp"

#include <stdexcept>
#include <string>

namespace MHD {

class Logger {
public:
    Logger();
    ~Logger();

    void LogError(Exception const& exception);

    void LogError(std::exception const& exception);

private:
    Error m_code;
    std::string m_message;
};

}

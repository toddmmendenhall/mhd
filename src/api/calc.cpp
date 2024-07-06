#include "calc.hpp"
#include "error.hpp"
#include "logger.hpp"
#include "safe.hpp"

#include <memory>

namespace MHD {

Calc::Calc() {
    m_logger = std::make_unique<Logger>();
    m_profile = std::make_unique<Profile>();
}

Calc::~Calc() {}

Error Calc::Run() {
    safe(m_logger, [&]() {
        // do stuff here
        return Error::SUCCESS;
    });
    return Error::SUCCESS;
}

}

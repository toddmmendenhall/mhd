#include "calc.hpp"
#include "error.hpp"

#include <memory>

namespace MHD {

Calc::Calc() {
    m_profile = std::make_unique<Profile>();
}

Calc::~Calc() {}

// Error Calc::Run() {
//     safe(m_logger, [&]() {
//         // do stuff here
//         return Error::SUCCESS;
//     });
//     return Error::SUCCESS;
// }

}

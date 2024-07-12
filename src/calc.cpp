#include "calc.hpp"
#include "error.hpp"

#include <memory>

namespace MHD {

Calc::Calc() {
    m_profile = std::make_unique<Profile>();
}

Calc::~Calc() {}

}

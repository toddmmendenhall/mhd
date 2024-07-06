#pragma once

#include "error.hpp"
#include "logger.hpp"

#include <functional>
#include <memory>

namespace MHD {

Error safe(std::unique_ptr<Logger> const& logger, std::function<Error()> const& unsafe);

}

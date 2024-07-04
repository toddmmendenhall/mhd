#pragma once

#include "exception.hpp"
#include "logger.hpp"

#include <functional>
#include <memory>

Exception api_safe(std::unique_ptr<Logger> const& logger, std::function<Exception()> api_unsafe);

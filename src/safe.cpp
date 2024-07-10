#include "error.hpp"
#include "exception.hpp"
#include "logger.hpp"
#include "safe.hpp"

#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace MHD {

Error safe(std::unique_ptr<Logger> const& logger, std::function<Error()> const& unsafe) {
    try {
        return unsafe();
    } catch (Exception const& exception) {
        logger->LogError(exception);
    } catch (std::exception const& exception) {
        logger->LogError(exception);
    } catch (...) {
        std::cout << "UNHANDLED EXCEPTION";
    }
    return Error::SUCCESS;
}

}

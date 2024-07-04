#include "api_safe.hpp"
#include "exception.hpp"
#include "logger.hpp"

#include <functional>
#include <memory>
#include <stdexcept>

Exception api_safe(std::unique_ptr<Logger> const& logger, std::function<Exception()> api_unsafe) {
    try {
        api_unsafe();
    } catch (Exception const& e) {
        logger.Error(e);
    } catch (std::exception const& e) {
        logger.Error(e);
    } catch (...) {
        logger.Error(UNHANDLED_EXCEPTION);
        std::cout << "UNHANDLED EXCEPTION";
    }
    return SUCCESS;
}

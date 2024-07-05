#include "error.hpp"
#include "exception.hpp"
#include "logger.hpp"
#include "safe.hpp"

#include <functional>
#include <memory>
#include <stdexcept>

Error safe(std::unique_ptr<Logger> const& logger, std::function<Error()> const unsafe) {
    try {
        return unsafe();
    } catch (Exception const& e) {
        logger.LogError(e);
    } catch (std::exception const& e) {
        logger.LogError(e);
    } catch (...) {
        std::cout << "UNHANDLED EXCEPTION";
    }
    return SUCCESS;
}

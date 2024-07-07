#pragma once

#include "error.hpp"

#include <string>

namespace MHD {

/**
 * \brief blah blah balh
 */
class Exception {
public:
    Exception();
    ~Exception();
    
    Error const GetError() const;
    std::string const& GetMessage() const;
    
    void SetError(Error const error);
    void SetMessage(std::string const& message);

private:
    Error m_code;
    std::string m_message;
};

}

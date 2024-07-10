#pragma once

#include "error.hpp"
#include "grid.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

/**
 * \brief Adds two numbers.
 *
 * This function takes two numbers, adds them, and then returns the result.
 * 
 * \param x The first number to add.
 * \param y The second number to add.
 * \return The sum of the two numbers.
 */
class Calc {
public:
    Calc();
    ~Calc();

    Error Run();
        
private:
    std::unique_ptr<Profile> m_profile;
    std::unique_ptr<Grid> m_grid;
};

}

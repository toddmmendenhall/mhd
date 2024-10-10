#pragma once

#include <memory>

namespace MHD {

class Grid;
class Profile;

class Calc {
public:
    Calc();
    ~Calc();

    void Run();
    
private:
    std::unique_ptr<Grid> m_grid;
    std::unique_ptr<Profile> m_profile;
};

} // namespace MHD

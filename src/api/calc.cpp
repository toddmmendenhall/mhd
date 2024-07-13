#include "calc.hpp"
#include "grid_factory.hpp"
#include "profile.hpp"

#include <iostream>

namespace MHD {

Calc::Calc() {
    Profile* m_profile = new Profile();
    GridFactory* gridFactory = new GridFactory();
    m_grid = gridFactory->CreateGrid(m_profile);
    delete gridFactory;
}

Calc::~Calc() {
    delete m_profile;
    delete m_grid;
}

void Calc::Run() {
    std::cout << "Running calc..." << std::endl;
}

Profile* Calc::GetProfile() const {return m_profile;}
Grid* Calc::GetGrid() const {return m_grid;}

} // namespace MHD

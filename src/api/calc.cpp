#include <calc.hpp>
#include <constants.hpp>
#include <execution_controller.hpp>
#include <grid.hpp>
#include <profile.hpp>
#include <solver.hpp>
#include <variable_store.hpp>

#include <iostream>
#include <fstream>
#include <memory>
#include <string>

namespace MHD {

Calc::Calc(Profile const& profile) : m_profile(profile) {
    m_executionController = std::make_unique<ExecutionController>();
    m_variableStore = std::make_unique<VariableStore>();
    m_grid = gridFactory(m_profile);
    m_solver = solverFactory(m_profile, *m_executionController, *m_variableStore, *m_grid);
}

Calc::~Calc() = default;

void Calc::SetInitialConditions() {
    std::size_t const size = 10;
    m_variableStore->rho.resize(size, ATMOSPHERIC_DENSITY_STP);
    m_variableStore->rhoU.resize(size, 0.0);
    m_variableStore->rhoV.resize(size, 0.0);
    m_variableStore->rhoW.resize(size, 0.0);
    m_variableStore->rhoE.resize(size, 0.0);
    m_variableStore->u.resize(size, 0.0);
    m_variableStore->v.resize(size, 0.0);
    m_variableStore->w.resize(size, 0.0);
    m_variableStore->p.resize(size, ATMOSPHERIC_PRESSURE_STP);
    m_variableStore->e.resize(size, ATMOSPHERIC_PRESSURE_STP);
    m_variableStore->bX.resize(size, 0.0);
    m_variableStore->bY.resize(size, 0.0);
    m_variableStore->bZ.resize(size, 0.0);
    m_variableStore->rhoInv.resize(size,0.0);
    m_variableStore->uu.resize(size, 0.0);
    m_variableStore->t.resize(size, ATMOSPHERIC_STANDARD_TEMPERATURE);
    m_variableStore->bb.resize(size, 0.0);
    m_variableStore->faceBX.resize(size, 0.0);
    m_variableStore->faceBY.resize(size, 0.0);
    m_variableStore->faceBZ.resize(size, 0.0);
    m_variableStore->edgeEX.resize(size, 0.0);
    m_variableStore->edgeEY.resize(size, 0.0);
    m_variableStore->edgeEZ.resize(size, 0.0);
}

void Calc::SetSodShockTube() {
    std::size_t const size = m_variableStore->numCells;
    for (std::size_t i = 0; i < size; ++i) {
        if (i >= size / 2) {
            m_variableStore->rho[i] *= 0.125;
            m_variableStore->p[i] *= 0.1;
        }
    }
}

void Calc::Run() {
    while (m_currentTime < m_duration) {
        std::cout << m_currentTime;
        if (OutputDataOption::YES == m_profile.m_outputDataOption) {
            WriteData(*m_variableStore);
        }
        m_solver->PerformTimeStep();
        m_currentTime += 0.01;
        m_currentStep++;
    }
}

void Calc::WriteData(VariableStore const& varStore) {
    std::ofstream myFile;
    std::string filename = "results.csv";
    
    // Open the file for writing
    myFile.open(filename);

    // Check if the file was successfully opened
    if (myFile.is_open()) {
        // Write data to the file
        myFile << "# rho, u, v, w, p, T" << std::endl;

        for (std::size_t i = 0; i < varStore.numCells; ++i) {
            myFile <<
            varStore.rho[i] << ", " <<
            varStore.u[i] << ", " <<
            varStore.v[i] << ", " <<
            varStore.w[i] << ", " <<
            varStore.p[i] << ", " <<
            varStore.t[i] << std::endl;
        }

        // Close the file
        myFile.close();
        std::cout << "Data written to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

} // namespace MHD

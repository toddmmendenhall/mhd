#include <calc.hpp>
#include <constants.hpp>
#include <execution_controller.hpp>
#include <grid.hpp>
#include <integration.hpp>
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
    m_solver = solverFactory(m_profile, *m_executionController, *m_variableStore);
    m_integrator = integratorFactor();
}

Calc::~Calc() = default;

void Calc::SetInitialConditions() {
    std::size_t const size = 10;
    m_variableStore->m_rho.resize(size, ATMOSPHERIC_DENSITY_STP);
    m_variableStore->m_rhoU.resize(size, 0.0);
    m_variableStore->m_rhoV.resize(size, 0.0);
    m_variableStore->m_rhoW.resize(size, 0.0);
    m_variableStore->m_rhoE.resize(size, 0.0);
    m_variableStore->m_u.resize(size, 0.0);
    m_variableStore->m_v.resize(size, 0.0);
    m_variableStore->m_w.resize(size, 0.0);
    m_variableStore->m_p.resize(size, ATMOSPHERIC_PRESSURE_STP);
    m_variableStore->m_e.resize(size, ATMOSPHERIC_PRESSURE_STP);
    m_variableStore->m_bX.resize(size, 0.0);
    m_variableStore->m_bY.resize(size, 0.0);
    m_variableStore->m_bZ.resize(size, 0.0);
    m_variableStore->m_rhoInv.resize(size,0.0);
    m_variableStore->m_uu.resize(size, 0.0);
    m_variableStore->m_t.resize(size, ATMOSPHERIC_STANDARD_TEMPERATURE);
    m_variableStore->m_bb.resize(size, 0.0);
    m_variableStore->m_faceBX.resize(size, 0.0);
    m_variableStore->m_faceBY.resize(size, 0.0);
    m_variableStore->m_faceBZ.resize(size, 0.0);
    m_variableStore->m_edgeEX.resize(size, 0.0);
    m_variableStore->m_edgeEY.resize(size, 0.0);
    m_variableStore->m_edgeEZ.resize(size, 0.0);
}

void Calc::SetSodShockTube() {
    std::size_t const size = m_variableStore->m_rho.size();
    for (std::size_t i = 0; i < size; ++i) {
        if (i >= size / 2) {
            m_variableStore->m_rho[i] *= 0.125;
            m_variableStore->m_p[i] *= 0.1;
        }
    }
}

void Calc::Run() {
    while (m_currentTime < m_duration) {
        std::cout << m_currentTime;
        if (OutputDataOption::YES == m_profile.m_outputDataOption) {
            WriteData(*m_variableStore);
        }
        m_solver->PerformTimeStep(*m_grid);
        m_integrator->Solve(*m_executionController, *m_variableStore);
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

        for (std::size_t i = 0; i < varStore.m_rho.size(); ++i) {
            myFile <<
            varStore.m_rho[i] << ", " <<
            varStore.m_u[i] << ", " <<
            varStore.m_v[i] << ", " <<
            varStore.m_w[i] << ", " <<
            varStore.m_p[i] << ", " <<
            varStore.m_t[i] << std::endl;
        }

        // Close the file
        myFile.close();
        std::cout << "Data written to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

} // namespace MHD

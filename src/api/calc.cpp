#include <calc.hpp>
#include <constants.hpp>
#include <error.hpp>
#include <execution_controller.hpp>
#include <grid.hpp>
#include <profile.hpp>
#include <solver.hpp>
#include <variable_store.hpp>

#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <string>

namespace MHD {

Calc::Calc(Profile const& profile) : m_profile(profile) {
    m_executionController = std::make_unique<ExecutionController>();
    m_grid = gridFactory(m_profile);
    m_variableStore = std::make_unique<VariableStore>(*m_grid);
    m_solver = solverFactory(m_profile, *m_executionController, *m_variableStore, *m_grid);
}

Calc::~Calc() = default;

void Calc::SetInitialCondition(InitialCondition ic) {
    if (InitialCondition::ATMOSPHERE == ic) {
        SetAtmosphere();
    } else if (InitialCondition::SOD_SHOCK_TUBE == ic) {
        SetSodShockTube();
    } else {
        throw Error::INVALID_INITIAL_CONDITION;
    }
}

void Calc::SetAtmosphere() {
    double const p = STANDARD_PRESSURE;
    double const rho = ATMOSPHERIC_DENSITY_STP;
    double const gamma = 1.4;
    double const e = p / ((gamma - 1.0) * rho);

    for (std::size_t i = 0; i < m_grid->NumCells(); ++i) {
        m_variableStore->rho[i] = rho;
        m_variableStore->rhoU[i] = 0.0;
        m_variableStore->rhoV[i] = 0.0;
        m_variableStore->rhoW[i] = 0.0;
        m_variableStore->rhoE[i] = rho * e;
    }
}

void Calc::SetSodShockTube() {
    double const gamma = 1.4;

    double const p1 = STANDARD_PRESSURE;
    double const rho1 = ATMOSPHERIC_DENSITY_STP;
    double const e1 = p1 / ((gamma - 1.0) * rho1);

    double const p2 = 0.1 * STANDARD_PRESSURE;
    double const rho2 = 0.125 * ATMOSPHERIC_DENSITY_STP;
    double const e2 = p2 / ((gamma - 1.0) * rho2);

    std::size_t const numCells = m_grid->NumCells();
    for (std::size_t i = 0; i < numCells; ++i) {
        if (i <= numCells / 2) {
            m_variableStore->rho[i] = rho1;
            m_variableStore->rhoU[i] = 0.0;
            m_variableStore->rhoV[i] = 0.0;
            m_variableStore->rhoW[i] = 0.0;
            m_variableStore->rhoE[i] = rho1 * e1;
        } else {
            m_variableStore->rho[i] = rho2;
            m_variableStore->rhoU[i] = 0.0;
            m_variableStore->rhoV[i] = 0.0;
            m_variableStore->rhoW[i] = 0.0;
            m_variableStore->rhoE[i] = rho2 * e2;
        }
    }
}

void Calc::Run() {
    while (m_currentTime < m_duration) {
        if (OutputDataOption::YES == m_profile.m_outputDataOption) {
            if (m_currentStep % 10 == 0) {
                WriteData(*m_variableStore);
                ++m_currentOutput;
            }
        }
        m_solver->PerformTimeStep();
        m_currentTime += m_solver->TimeStep();
        m_currentStep++;
    }
}

void Calc::WriteData(VariableStore const& varStore) {
    std::ofstream myFile;
    std::string filename = "results_" + std::to_string(m_currentOutput) + ".csv";
    
    // Open the file for writing
    myFile.open(filename);

    // Check if the file was successfully opened
    if (myFile.is_open()) {
        // Write data to the file
        myFile << "# x, rho, u, v, w, p, e, cc, T" << std::endl;
        myFile << "# time: " << m_currentTime << " s" << std::endl;

        for (std::size_t i = 0; i < m_grid->NumNodes(); ++i) {
            myFile <<
            m_grid->Nodes()[i][0] << ", " <<
            varStore.rho[i] << ", " <<
            varStore.u[i] << ", " <<
            varStore.v[i] << ", " <<
            varStore.w[i] << ", " <<
            varStore.p[i] << ", " <<
            varStore.e[i] << ", " <<
            varStore.cs[i] << ", " <<
            varStore.t[i] << std::endl;
        }

        // Close the file
        myFile.close();
        std::cout << "Data written to " << filename.data() << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

} // namespace MHD

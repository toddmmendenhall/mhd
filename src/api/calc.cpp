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
#include <memory>
#include <string>

namespace MHD {

Calc::Calc(Profile const& profile) : m_profile(profile) {
    m_executionController = std::make_unique<ExecutionController>();
    m_grid = gridFactory(m_profile);
    m_variableStore = std::make_unique<VariableStore>(*m_grid);
    m_solver = solverFactory(m_profile, *m_executionController, *m_variableStore, *m_grid, tStep);
}

Calc::~Calc() = default;

void Calc::SetInitialCondition(InitialCondition ic) {
    if (InitialCondition::ATMOSPHERE == ic) {
        return;
    } else if (InitialCondition::SOD_SHOCK_TUBE == ic) {
        SetSodShockTube();
    } else {
        throw Error::INVALID_INITIAL_CONDITION;
    }
}

void Calc::SetSodShockTube() {
    std::size_t const size = m_grid->NumInteriorCells();
    for (std::size_t i = 0; i < size; ++i) {
        if (i >= size / 2) {
            m_variableStore->rho[i] *= 0.125;
            m_variableStore->p[i] *= 0.1;
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
        m_currentTime += tStep;
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
        myFile << "# rho, u, v, w, p, T" << std::endl;

        for (std::size_t i = 0; i < m_grid->NumCells(); ++i) {
            myFile <<
            varStore.rho[i] << ", " <<
            varStore.u[i] << ", " <<
            varStore.v[i] << ", " <<
            varStore.w[i] << ", " <<
            varStore.p[i] << ", " <<
            varStore.e[i] << ", " <<
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

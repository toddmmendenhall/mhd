#include <calc.hpp>
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
    m_solver = solverFactory(m_profile, *m_executionController, *m_variableStore);
}

Calc::~Calc() = default;

void Calc::Run() {
    if (OutputDataOption::YES == m_profile.m_outputDataOption) {
        WriteData(*m_variableStore);
    }
    m_solver->PerformTimeStep();
}

void Calc::WriteData(VariableStore const& varStore) {
    std::ofstream myFile;
    std::string filename = "example.txt";
    
    // Open the file for writing
    myFile.open(filename);

    // Check if the file was successfully opened
    if (myFile.is_open()) {
        // Write data to the file
        myFile << "Hello, file!" << std::endl;
        myFile << "This is another line." << std::endl;

        // Close the file
        myFile.close();
        std::cout << "Data written to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

} // namespace MHD

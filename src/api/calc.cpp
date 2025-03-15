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
    } else if (InitialCondition::BRIO_WU_SHOCK_TUBE == ic) {
        SetBrioWuShockTube();
    } else {
        throw Error::INVALID_INITIAL_CONDITION;
    }
}

void Calc::SetAtmosphere() {
    double const rho = ATMOSPHERIC_DENSITY_STP;
    double const p = STANDARD_PRESSURE;
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

    double const rho1 = 1.0;
    double const p1 = 1.0 * STANDARD_PRESSURE;
    double const e1 = p1 / ((gamma - 1.0) * rho1);

    double const rho2 = 0.125;
    double const p2 = 0.1 * STANDARD_PRESSURE;
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

void Calc::SetBrioWuShockTube() {
    double const gamma = 2.0;
    double const bx = 0.75;
    double const bz = 0.0;

    double const rho1 = 1.0;
    double const p1 = 1.0 * STANDARD_PRESSURE;
    double const by1 = 1.0;
    double const b1Squared = bx * bx + by1 * by1 + bz * bz;
    double const e1 = p1 / ((gamma - 1.0) * rho1);

    double const rho2 = 0.125;
    double const p2 = 0.1 * STANDARD_PRESSURE;
    double const by2 = -1.0;
    double const b2Squared = bx * bx + by2 * by2 + bz * bz;
    double const e2 = p2 / ((gamma - 1.0) * rho2);

    std::size_t const numCells = m_grid->NumCells();
    for (std::size_t i = 0; i < numCells; ++i) {
        if (i <= numCells / 2) {
            m_variableStore->rho[i] = rho1;
            m_variableStore->rhoU[i] = 0.0;
            m_variableStore->rhoV[i] = 0.0;
            m_variableStore->rhoW[i] = 0.0;
            m_variableStore->rhoE[i] = rho1 * e1 + 0.5 * b1Squared;
            m_variableStore->bx[i] = bx;
            m_variableStore->by[i] = by1;
            m_variableStore->bz[i] = bz;
        } else {
            m_variableStore->rho[i] = rho2;
            m_variableStore->rhoU[i] = 0.0;
            m_variableStore->rhoV[i] = 0.0;
            m_variableStore->rhoW[i] = 0.0;
            m_variableStore->rhoE[i] = rho2 * e2 + 0.5 * b2Squared;
            m_variableStore->bx[i] = bx;
            m_variableStore->by[i] = by2;
            m_variableStore->bz[i] = bz;
        }
    }
}

void Calc::Run() {
    while (m_currentTime < m_duration) {
        m_solver->PrimFromCons();
        if (OutputDataOption::YES == m_profile.m_outputDataOption) {
            if (m_currentTime >= m_currentOutput * m_outputPeriod) {
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
        myFile << "# x, rho, u, v, w, bx, by, bz, e, p, T, cs" << std::endl;
        myFile << "# time: " << m_currentTime << " s" << std::endl;

        for (std::size_t i = 0; i < m_grid->NumNodes(); ++i) {
            myFile <<
            m_grid->Nodes()[i][0] << ", " <<
            varStore.rho[i] << ", " <<
            varStore.u[i] << ", " <<
            varStore.v[i] << ", " <<
            varStore.w[i] << ", " <<
            varStore.bx[i] << ", " <<
            varStore.by[i] << ", " <<
            varStore.bz[i] << ", " <<
            varStore.e[i] << ", " <<
            varStore.p[i] << ", " <<
            varStore.t[i] << ", " <<
            varStore.cs[i] << std::endl;
        }

        // Close the file
        myFile.close();
        std::cout << "Data written at time: " << m_currentTime << " s" << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

} // namespace MHD

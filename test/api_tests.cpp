#include <calc.hpp>
#include <profile.hpp>
#include <profile_options.hpp>

#include "gtest/gtest.h"

#include <iostream>

TEST(APITests, CanCreateProfile) {
    MHD::Profile profile;
    EXPECT_EQ(MHD::Dimension::ONE, profile.m_gridDimensionOption);
}

TEST(APITests, CanCreateCalc) {
    MHD::Profile profile;
    profile.m_outputDataOption = MHD::OutputDataOption::YES;
    MHD::Calc calc(profile);
    calc.SetInitialConditions();
    calc.Run();
}

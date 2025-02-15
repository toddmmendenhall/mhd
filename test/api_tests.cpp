#include <profile.hpp>
#include <profile_options.hpp>

#include "gtest/gtest.h"

#include <iostream>

TEST(APITests, CanCreateProfile) {
    MHD::Profile profile;
    EXPECT_EQ(MHD::Dimension::ONE, profile.m_gridDimensionOption);
}

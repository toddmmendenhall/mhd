#pragma once

namespace MHD {

enum class Dimension {
    ONE = 0,
    TWO = 1,
    THREE = 2,
};

enum class Geometry {
    CARTESIAN = 0,
};

enum class BoundaryCondition {
    NOT_A_BOUNDARY = 0,
    DIRICHLET = 1,
    NEUMANN = 2,
};

enum class SpatialDerivativeMethod {
    FINITE_DIFFERENCE = 0,
    FINITE_VOLUME = 1,
};

enum class TemporalIntegrationMethod {
    FORWARD_EULER = 0,
};

enum class EquationOfState {
    CALORICALLY_PERFECT_GAS = 0,
    THERMALLY_PERFECT_GAS = 1,
};

enum class FluxScheme {
    HLLC = 0,
    KT = 1,
    HOG = 2,    // High Order Goduov
    LOG = 3,    // Low Order Godunov
};

enum class ReconstructionOption {
    LINEAR = 0,
};

} // namespace MHD

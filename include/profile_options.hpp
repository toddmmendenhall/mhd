#pragma once

namespace MHD {

enum class Dimension {
    ONE = 0,
    TWO = 1,
    THREE = 2,
};

enum class BoundaryConditionOption {
    REFLECTIVE = 0,
    OUTFLOW = 1,
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
    CONSTANT = 0,
    LINEAR = 1,
};

enum class CompressibleOption {
    COMPRESSIBLE = 0,
    INCOMPRESSIBLE = 1,
};

enum class ViscousOption {
    VISCOUS = 0,
    INVISCID = 1,
};

enum class HydroOption {
    MAGNETO_HYDRO = 0,
    PURE_HYDRO = 1,
};

enum class OutputDataOption {
    NO = 0,
    YES = 1,
};

} // namespace MHD

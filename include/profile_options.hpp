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

enum class TemporalIntegrationMethod {
    FORWARD_EULER = 0,
};

enum class EquationOfState {
    CALORICALLY_PERFECT_GAS = 0,
    THERMALLY_PERFECT_GAS = 1,
};

enum class FluxScheme {
    KT = 0,
};

enum class ReconstructionOption {
    CONSTANT = 0,
    LINEAR = 1,
};

enum class CompressibleOption {
    COMPRESSIBLE = 0,
};

enum class OutputDataOption {
    NO = 0,
    YES = 1,
};

} // namespace MHD

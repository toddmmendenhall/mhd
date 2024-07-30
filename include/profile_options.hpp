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

} // namespace MHD

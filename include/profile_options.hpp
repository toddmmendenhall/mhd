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
    DIRICHLET = 0,
    NEUMANN = 1,
};

} // namespace MHD
#pragma once

namespace MHD {

enum class Error {
    SUCCESS = 0,
    INVALID_FLUX_SCHEME = 1,
    INVALID_GRID_DIMENSION = 2,
    INVALID_GRID_GEOMETRY = 3,
    INVALID_RECONSTRUCTION_OPTION = 4,
};

}

# Build instructions
* Clone this repository:
* Update submodules: `git update submodules --init`
* Create a build directory: `mkdir build`
* Change directories to the build directory: `cd build`
* Run the cmake generate step: `cmake .. -DCMAKE_BUILD_TYPE=Debug`
* Run the cmake build step: `cmake --build . --config Debug`
* Run tests: `./test/mhd_tests`

# Main user-facing objects
* Profile - Contains all of the options for running a simulation
* Calc - Top level calculation object for simulation control flow

# Main private objects owned by Calc
* Grid - Polymorphic class depending on the user-specified options
* Solver - Polymorphic class depending on the user-specified options
* VariableStore - Holds all the data necessary to run a calc
* ExecutionController - Responsible for launching and cleaning up kernels

# Main private objects own by Solver
* IReconstruction - Polymorphic class responsible for computing the left and right states of each cell
* IFlux - Polymorphic class responsible for computing the fluxes
* ITransport - Polymorphic class responsbile for computing the residuals for transport terms
* ISource - Polymorphic class responsible for computing the residuals for source terms
* IIntegrator - Polymorphic class responsible for computing the time integration

# Naming conventions
* Filenames - snake_case
* Namespaces - ALL_CAPS_WITH_UNDERSCORES
* Classes - UpperCamelCase
* Methods - UpperCamelCase
* Free functions - lowerCamelCase
* Member classes/structs - m_lowerCamelCase
* Member plain old data variables - lowerCamelCase
* Arguments - lowerCamelCase
* Constants - ALL_CAPS_WITH_UNDERSCORES
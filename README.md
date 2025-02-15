# Main user-facing objects
* Profile - Contains all of the options for running a simulation
* Calc - Top level calculation object for simulation control flow

# Main private objects owned by Calc
* Grid - Polymorphic class depending on the user-specified options
* Solver - Polymorphic class depending on the user-specified options
* VariableStore - Holds all the data necessary to run a calc
* ExecutionController - Responsible for launching and cleaning up kernels

# Naming conventions
* Filenames - snake_case
* Namespaces - ALL_CAPS_WITH_UNDERSCORES
* Classes - UpperCamelCase
* Methods - UpperCamelCase
* Free functions - lowerCamelCase
* Member variables - m_lowerCamelCase
* Arguments - lowerCamelCase
* Constants - ALL_CAPS_WITH_UNDERSCORES
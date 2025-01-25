#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
class Profile;
class VariableStore;

class Solver {
public:
    Solver(Profile const& profile);
    ~Solver();

    void ComputePrimitivesFromConserved();

    void UpdateConservedFromPrimitives();

private:
    std::unique_ptr<VariableStore> varStore;
    std::unique_ptr<ExecutionController> execCtrl;
};

} // namespace MHD

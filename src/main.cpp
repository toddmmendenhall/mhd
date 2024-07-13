#include "calc.hpp"
#include "profile.hpp"

int main() {
    MHD::Calc* calc = new MHD::Calc();
    calc->Run();
    delete calc;
    return 0;
}

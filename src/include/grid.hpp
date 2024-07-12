#pragma once

namespace MHD {

class Grid {
public:
    Grid();
    ~Grid();

    virtual void Method() = 0;
};

}

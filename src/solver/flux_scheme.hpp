#pragma once

#include <memory>

namespace MHD {

class Face;
class Profile;

class IFluxScheme {
public:
    IFluxScheme() = default;
    virtual ~IFluxScheme() = default;
    virtual void blah() const = 0;
    virtual void operator()(std::size_t const idx) {}

};

std::unique_ptr<IFluxScheme> flux_scheme_factory(Profile const& profile, std::size_t const numFaces);

} // namespace MHD
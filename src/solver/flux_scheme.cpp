#include <error.hpp>
#include <flux_scheme.hpp>
#include <profile.hpp>
#include <profile_options.hpp>
#include <face.hpp>

namespace MHD {

struct HLLCKernel : IFluxScheme {
    HLLCKernel(std::size_t const numFaces) {}
    void blah() const {}

    inline void operator()(std::size_t const idx) {}
};

class HLLC : public IFluxScheme {
public:
    HLLC(std::size_t const numFaces) : IFluxScheme() {}
    void blah() const {}
};

std::unique_ptr<IFluxScheme> flux_scheme_factory(Profile const& profile, std::size_t const numFaces) {
    auto fluxSchemeName = profile.GetFluxScheme();
    switch (fluxSchemeName) {
    case FluxScheme::HLLC:
        return std::make_unique<HLLCKernel>(numFaces);
    default:
        throw Error::INVALID_FLUX_SCHEME;
    }
}

} // namespace MHD
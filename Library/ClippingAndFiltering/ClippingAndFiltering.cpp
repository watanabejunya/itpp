#include <itpp/itbase.h>
#include "../../Library/Common/helper.hpp"
#include "../../Library/ClippingAndFiltering/ClippingAndFiltering.hpp"

using std::acos;
using std::exp;
using std::pow;
using std::erfc;
using itpp::cvec;
using helper::get_average_power;
using helper::oversample;
using helper::downsample;
using helper::fft_normalized;
using helper::ifft_normalized;

const double PI = acos(-1);

ClippingAndFiltering::ClippingAndFiltering(double clippingRatio, int oversamplingFactor)
: clippingRatio(clippingRatio), oversamplingFactor(oversamplingFactor)
{
}

ClippingAndFiltering::~ClippingAndFiltering()
{
}

double ClippingAndFiltering::calculate_attenuation_factor()
{
    switch (this->clipModel) {
        case CLIP_MODEL::SoftLimiter:
        return 1.0 - exp(-pow(this->clippingRatio, 2.0)) + sqrt(PI) * this->clippingRatio * erfc(this->clippingRatio) / 2.0;
    default:
        return 1.0 - exp(-pow(this->clippingRatio, 2.0)) + sqrt(PI) * this->clippingRatio * erfc(this->clippingRatio) / 2.0;
    };
}

cvec ClippingAndFiltering::clip(const cvec& signals)
{
    switch (this->clipModel) {
        case CLIP_MODEL::SoftLimiter:
        return clip_by_soft_limitter(signals);
    default:
        return clip_by_soft_limitter(signals);
    };
}

cvec ClippingAndFiltering::clip_by_soft_limitter(const cvec& signals)
{
    cvec clipped = signals;
    double maxAmplitude;

    // Configure threshold amplitude.
    maxAmplitude = this->clippingRatio * sqrt(get_average_power(signals));

    for (int i = 0, size = signals.size(); i < size; i++) {
        if (abs(clipped[i]) > maxAmplitude) {
            clipped[i] *= maxAmplitude / abs(clipped[i]);
        }
    }

    return clipped;
}

cvec ClippingAndFiltering::normalize_attenuation(const cvec& signals)
{
    double attenuationFactor = this->calculate_attenuation_factor();

    return signals / attenuationFactor;
}

cvec ClippingAndFiltering::clip_and_filter(const cvec& symbols)
{
    cvec signals = ifft_normalized(oversample(symbols, this->oversamplingFactor));

    signals = this->clip(signals);

    return downsample(fft_normalized(signals), this->oversamplingFactor);

}

#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include <itpp/signal/transforms.h>
#include "../../Library/Channel/StaticFadingChannel.hpp"

using itpp::linspace_fixed_step;
using itpp::exp;
using itpp::sqrt;
using itpp::norm;
using itpp::vec;
using itpp::ivec;
using itpp::cvec;
using itpp::ones;
using itpp::zeros_c;
using itpp::elem_mult;
using itpp::fft;


StaticFadingChannel::StaticFadingChannel(int numberOfTaps, const DelayProfileType& type)
{
    this->numberOfTaps = numberOfTaps;
    this->delayProfile = linspace_fixed_step(0, numberOfTaps - 1);
    this->impulseResponses.set_size(numberOfTaps);

    switch (type) {
        case DelayProfileType::Uniform:
            this->amplitudeProfile = sqrt(ones(numberOfTaps) / (double)numberOfTaps);
            this->amplitudeProfile /= norm(this->amplitudeProfile);
        case DelayProfileType::Exponential:
            this->amplitudeProfile = sqrt(exp(linspace_fixed_step(0.0, - numberOfTaps + 1.0, -1.0)));
            this->amplitudeProfile /= norm(this->amplitudeProfile);
    }
}

StaticFadingChannel::~StaticFadingChannel()
{
}

void StaticFadingChannel::set_channel_profile(const vec& amplitudeProfile, const ivec& delayProfile)
{
    it_assert_debug(min(delayProfile) == 0,
             "StaticFadingChannel::set_channel_profile(): Minimum relative delay must be 0.");
    it_assert_debug(amplitudeProfile.size() == delayProfile.size(),
             "StaticFadingChannel::set_channel_profile(): Power and delay vectors must be of equal length!");
    it_assert_debug(delayProfile(0) == 0,
             "StaticFadingChannel::set_channel_profile(): First tap must be at zero delay");
    for (int i = 1; i < delayProfile.size(); i++) {
        it_assert_debug(delayProfile[i] > delayProfile[i - 1],
            "StaticFadingChannel::set_channel_profile(): Delays should be sorted and unique");
    }

    this->numberOfTaps = delayProfile.size();
    this->amplitudeProfile = amplitudeProfile;
    this->delayProfile = delayProfile;
}

void StaticFadingChannel::set_fft_size(int fftSize, int oversamplingFactor)
{
    this->fftSize = fftSize;
    this->oversamplingFactor = oversamplingFactor;
}

void StaticFadingChannel::init()
{
    for (int i = 0; i < this->numberOfTaps; i++) {
        this->impulseResponses[i] = this->generator.sample() * this->amplitudeProfile[i];
    }
}

cvec StaticFadingChannel::operator()(const cvec& input)
{
    int maxDelay = this->delayProfile[this->numberOfTaps - 1];
    cvec output = zeros_c(input.size() + maxDelay);

    for (int i = 0; i < this->numberOfTaps; i++) {
        output += concat(
            zeros_c(this->delayProfile[i]),
            input * this->impulseResponses[i],
            zeros_c(maxDelay - this->delayProfile[i])
        );
    }

    return output.left(input.size());
}

cvec StaticFadingChannel::get_channel_cofficients()
{
    it_assert(fftSize > this->delayProfile[this->numberOfTaps - 1],
        "StaticFadingChannel::get_channel_cofficients(): FFT size must be larger than the maximum delay.");
    it_assert(fftSize > this->delayProfile[this->numberOfTaps - 1],
        "StaticFadingChannel::get_channel_cofficients(): oversampling factor must be >= 0.");

    cvec impulseResponses = concat(
        this->impulseResponses,
        zeros_c(this->oversamplingFactor * this->fftSize - this->numberOfTaps)
    );

    cvec channelCofficients = fft(impulseResponses);

    return concat(
        channelCofficients.right(this->fftSize / 2),
        channelCofficients.left(this->fftSize / 2)
    );
}

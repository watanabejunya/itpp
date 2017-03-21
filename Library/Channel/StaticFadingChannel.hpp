#ifndef StaticFadingChannel_hpp
#define StaticFadingChannel_hpp

#include <itpp/itbase.h>
#include <itpp/base/random.h>

enum class DelayProfileType {Uniform = 1, Exponential = 2};

class StaticFadingChannel
{
    private:
        int numberOfTaps;
        int fftSize = 0;
        int oversamplingFactor = 1;
        itpp::vec amplitudeProfile;
        itpp::ivec delayProfile;
        itpp::cvec impulseResponses;
        itpp::cvec channelCoefficients;
        itpp::Complex_Normal_RNG generator;

        void set_channel_profile(const itpp::vec&, const itpp::ivec&);

    public:
        StaticFadingChannel(int, const DelayProfileType&);
        ~StaticFadingChannel();
        void set_fft_size(int, int = 1);
        void init();
        itpp::cvec operator()(const itpp::cvec&);
        itpp::cvec get_channel_cofficients();

};

#endif

#ifndef OFDM_hpp
#define OFDM_hpp

#include <itpp/itbase.h>

class OFDM
{
    private:
        int numberOfSubcarriers;
        int guardInterval;
        int oversamplingFactor;
        double normFactor;

        itpp::cvec oversample(const itpp::cvec&);
        itpp::cvec downsample(const itpp::cvec&);

    public:
        OFDM(int, int = 0, int = 1);
        ~OFDM();
        int get_number_of_subcarriers() const;
        itpp::cvec modulate(const itpp::cvec&);
        itpp::cvec demodulate(const itpp::cvec&);
};

#endif

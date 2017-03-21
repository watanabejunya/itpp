#include <itpp/itbase.h>
#include <itpp/signal/transforms.h>
#include "../../Library/OFDM/OFDM.hpp"
#include "../../Library/Common/helper.hpp"

using std::sqrt;
using itpp::vec;
using itpp::zeros_c;
using itpp::concat;
using itpp::cvec;
using itpp::fft;
using itpp::ifft;
using helper::is_power_of_2;


OFDM::OFDM(int numberOfSubcarriers, int guardInterval, int oversamplingFactor)
{
    it_assert_debug(numberOfSubcarriers >= 2,
        "OFDM::ODM(): Number of subcarriers must be >= 2.");
    it_assert_debug(is_power_of_2(numberOfSubcarriers),
        "OFDM::ODM(): Number of subcarriers must be power of 2.");
    it_assert_debug(1 <= oversamplingFactor && oversamplingFactor <= 8,
            "OFDM::ODM(): Oversampling factor must be >= 1 and <= 8.");
    it_assert_debug(oversamplingFactor == 1 || is_power_of_2(oversamplingFactor),
        "OFDM::ODM(): Oversampling factor must be 1 or power of 2.");
    it_assert_debug(0 <= guardInterval && guardInterval <= numberOfSubcarriers,
        "OFDM::ODM(): Guard interval must be >= 0 and <= # of subcarriers.");

    this->numberOfSubcarriers = numberOfSubcarriers;
    this->oversamplingFactor = oversamplingFactor;
    this->guardInterval = guardInterval;
    this->normFactor = sqrt(oversamplingFactor * numberOfSubcarriers);
}

OFDM::~OFDM()
{
}

cvec OFDM::oversample(const cvec& input)
{
    return concat(
        input.right(this->numberOfSubcarriers / 2),
        zeros_c(this->numberOfSubcarriers * (this->oversamplingFactor - 1)),
        input.left(this->numberOfSubcarriers / 2)
    );
}

cvec OFDM::downsample(const cvec& input)
{
    return concat(
        input.right(this->numberOfSubcarriers / 2),
        input.left(this->numberOfSubcarriers / 2)
    );
}

int OFDM::get_number_of_subcarriers() const
{
    return this->numberOfSubcarriers;
}

itpp::cvec OFDM::modulate(const itpp::cvec& input)
{
    it_assert_debug(input.size() == numberOfSubcarriers,
        "OFDM::modulate(): Signal size does not match subbcarries.");

    cvec signal = ifft(this->oversample(input)) * this->normFactor;

    return concat(signal.right(this->oversamplingFactor * this->guardInterval), signal);
}

itpp::cvec OFDM::demodulate(const itpp::cvec& input)
{
    it_assert_debug(input.size() == this->oversamplingFactor * (this->numberOfSubcarriers + this->guardInterval),
        "OFDM::demodulate(): Signal size does not match subbcarries.");

    cvec signal = fft(input.right(this->oversamplingFactor * this->numberOfSubcarriers));

    return this->downsample(signal / this->normFactor);
}

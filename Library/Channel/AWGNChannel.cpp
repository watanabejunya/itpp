#include <itpp/itbase.h>
#include "../../Library/Channel/AWGNChannel.hpp"

using std::sqrt;
using itpp::sqr;
using itpp::cvec;


AWGNChannel::AWGNChannel()
{
}

AWGNChannel::~AWGNChannel()
{
}

void AWGNChannel::set_noise(double N0)
{
    it_assert_debug(N0 > 0, "AWGNChannel::set_noise(): Minimum relative delay must be 0.");

    this->sigma = sqrt(N0);
}

double AWGNChannel::get_noise() const
{
    return sqr(this->sigma);
}

cvec AWGNChannel::operator()(const cvec& input)
{
    cvec noise(input.size());
    cvec output(input.size());

    this->generator.sample_vector(input.size(), noise);
    noise *= sigma;

    output = input + noise;

    return output;
}

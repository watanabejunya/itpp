#include <itpp/itbase.h>
#include <itpp/comm/modulator.h>
#include "../../Library/Counter/SERCounter.hpp"
#include "../../Library/Common/helper.hpp"

using itpp::bvec;
using itpp::cvec;
using Modulator = itpp::Modulator<std::complex<double>>;
using helper::hard_decide_LLRs;


SERCounter::SERCounter(const Modulator& modulator)
{
    this->set_modulator(modulator);
}

SERCounter::~SERCounter()
{
}

void SERCounter::set_modulator(const Modulator& modulator)
{
    this->modulator = modulator;
    this->blockSize = modulator.get_k();
}

void SERCounter::count(const itpp::cvec& symbols1, const itpp::cvec& symbols2, double noise)
{
    it_assert_debug(symbols1.size() == symbols2.size(),
        "SERCounter::count(): The size of 2 symbol sequence must be the same.");

    bool isError;
    bvec transmittedBits = this->modulator.demodulate_bits(symbols1);
    bvec receivedBits = hard_decide_LLRs(this->modulator.demodulate_soft_bits(symbols2, noise));

    for (int i = 0; i < symbols1.size(); i++) {
        isError = false;

        for (int j = 0; j < this->blockSize; j++) {
            if (transmittedBits[this->blockSize * i + j] != receivedBits[this->blockSize * i + j]) {
                isError = true;
                break;
            }
        }

        if (isError) {
            this->errors++;
        } else {
            this->corrects++;
        }
    }
}

void SERCounter::count(const itpp::cvec& symbols1, const itpp::cvec& symbols2, const cvec& channel, double noise)
{
    it_assert_debug(symbols1.size() == symbols2.size(),
        "SERCounter::count(): The size of 2 symbol sequence must be the same.");

    bool isError;
    bvec transmittedBits = this->modulator.demodulate_bits(symbols1);
    bvec receivedBits = hard_decide_LLRs(this->modulator.demodulate_soft_bits(symbols2, channel, noise));

    for (int i = 0; i < symbols1.size(); i++) {
        isError = false;

        for (int j = 0; j < this->blockSize; j++) {
            if (transmittedBits[this->blockSize * i + j] != receivedBits[this->blockSize * i + j]) {
                isError = true;
                break;
            }
        }

        if (isError) {
            this->errors++;
        } else {
            this->corrects++;
        }
    }
}

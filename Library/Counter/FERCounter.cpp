#include <itpp/itbase.h>
#include <itpp/comm/modulator.h>
#include "../../Library/Counter/FERCounter.hpp"
#include "../../Library/Common/helper.hpp"

using itpp::bvec;
using itpp::cvec;
using Modulator = itpp::Modulator<std::complex<double>>;
using helper::hard_decide_LLRs;


FERCounter::FERCounter(int blockSize)
{
    it_assert_debug(blockSize > 0,
        "FERCounter::FERCounter(): The block size must be > 0.");

    this->blockSize = blockSize;
}

FERCounter::~FERCounter()
{
}

void FERCounter::count(const itpp::bvec& bits1, const itpp::bvec& bits2)
{
    it_assert_debug(bits1.size() == this->blockSize,
        "FERCounter::count(): The size of symbol sequence must be bigger than block size.");
    it_assert_debug(bits1.size() == bits2.size(),
        "FERCounter::count(): The size of 2 symbol sequence must be the same.");

    bool isError = false;

    for (int i = 0; i < this->blockSize; i++) {
        if (bits1[i] != bits2[i]) {
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

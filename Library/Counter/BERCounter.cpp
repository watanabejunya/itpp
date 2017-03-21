#include <itpp/itbase.h>
#include <itpp/comm/modulator.h>
#include "../../Library/Counter/BERCounter.hpp"
#include "../../Library/Common/helper.hpp"

using itpp::bvec;
using helper::hard_decide_LLRs;


BERCounter::BERCounter()
{
}

BERCounter::~BERCounter()
{
}

void BERCounter::count(const itpp::bvec& bits1, const itpp::bvec& bits2)
{
    it_assert_debug(bits1.size() == bits2.size(),
        "BERCounter::count(): The size of 2 bit sequence must be the same.");

    for (int i = 0; i < bits1.size(); i++) {
        if (bits1[i] == bits2[i]) {
            this->corrects++;
        } else {
            this->errors++;
        }
    }
}

#include <itpp/itbase.h>
#include "../../Library/Interleaver/Interleaver.hpp"
#include "../../Library/Interleaver/SRandomInterleaver.hpp"

using std::swap;
using itpp::randu;
using itpp::bvec;
using itpp::ivec;


SRandomInterleaver::SRandomInterleaver(const int size, const int s)
{
    this->size = size;
    this->s = (s > 0) ? s : floor(sqrt(size / 2.0));

    this->make_new_interleaver();
}

SRandomInterleaver::~SRandomInterleaver()
{
}

void SRandomInterleaver::make_new_interleaver()
{
    ivec interleaver = sort_index(randu(this->size));

    for (int i = 1, trial = 0, k = 0; i < this->size; i++, trial++) {
        bool failed = false;

        for (int j = 1; j <= this->s && j <= i; j++) {
            if (abs(interleaver[i + k] - interleaver[i - j]) <= this->s) {
                if (i + k < this->size - 1) {
                    failed = true;
                    k++;
                    i--;
                } else {
                    interleaver.shift_right(interleaver.right(this->size - i));
                    i = 0;
                    k = 0;
                }
                break;
            }
        }
        if (! failed) {
            swap(interleaver[i], interleaver[i + k]);
            k = 0;
        }

        // Reset first random data set.
        if (trial % (100 * this->size) == 0) {
            i = 0;
            k = 0;
            trial = 0;
            interleaver = sort_index(randu(this->size));
        }
    }

    this->interleaver = interleaver;
}

ivec SRandomInterleaver::make_interleaver(const int size, const int s)
{
    SRandomInterleaver interleaver(size, s);
    return interleaver.get_interleaver();
}

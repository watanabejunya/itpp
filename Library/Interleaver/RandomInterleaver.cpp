#include <itpp/itbase.h>
#include "../../Library/Interleaver/Interleaver.hpp"
#include "../../Library/Interleaver/RandomInterleaver.hpp"

using std::swap;
using itpp::randu;
using itpp::bvec;
using itpp::ivec;
using itpp::vec;


RandomInterleaver::RandomInterleaver(int size)
{
    this->size = size;

    this->make_new_interleaver();
}

RandomInterleaver::~RandomInterleaver()
{
}

void RandomInterleaver::make_new_interleaver()
{
    this->interleaver = sort_index(randu(this->size));
}

ivec RandomInterleaver::make_interleaver(const int size)
{
    RandomInterleaver interleaver(size);
    return interleaver.get_interleaver();
}

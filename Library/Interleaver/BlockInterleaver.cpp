#include <itpp/itbase.h>
#include "../../Library/Interleaver/Interleaver.hpp"
#include "../../Library/Interleaver/BlockInterleaver.hpp"

using std::swap;
using itpp::randu;
using itpp::bvec;
using itpp::ivec;


BlockInterleaver::BlockInterleaver(const int rows, const int columns)
{
    this->rows = rows;
    this->columns = columns;
    this->size = rows * columns;

    this->make_new_interleaver();
}

BlockInterleaver::~BlockInterleaver()
{
}

void BlockInterleaver::make_new_interleaver()
{
    this->interleaver.set_size(this->size);

    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->columns; j++) {
            this->interleaver[j * this->rows + i] = i * this->columns + j;
        }
    }
}

ivec BlockInterleaver::make_interleaver(const int size, const int s)
{
    BlockInterleaver interleaver(size, s);
    return interleaver.get_interleaver();
}

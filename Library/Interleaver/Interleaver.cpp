#include <itpp/itbase.h>
#include "../../Library/Common/helper.hpp"
#include "../../Library/Interleaver/Interleaver.hpp"

using std::swap;
using itpp::randu;
using itpp::bvec;
using itpp::ivec;
using itpp::vec;

Interleaver::Interleaver()
{
}

Interleaver::~Interleaver()
{
}

void Interleaver::make_new_interleaver()
{
}

int Interleaver::get_size() const
{
    return this->size;
}

ivec Interleaver::get_interleaver() const
{
    return this->interleaver;
}

void Interleaver::set_interleaver(const itpp::ivec& interleaver)
{
    this->size = interleaver.size();
    this->interleaver = interleaver;
}

bvec Interleaver::interleave(const bvec& input)
{
    int inputSize = input.size();
    int numberOfBlock = inputSize / this->size;
    bvec output(inputSize);

    it_assert_debug(inputSize % this->size == 0, "Interleaver::interleave(): Invalid length of given vector.");

    for (int i = 0; i < numberOfBlock; i++) {
        for (int j = 0; j < this->size; j++) {
            output[i * this->size + j] = input[i * this->size + this->interleaver[j]];
        }
    }

    return output;
}

vec Interleaver::interleave(const vec& input)
{
    int inputSize = input.size();
    int numberOfBlock = inputSize / this->size;
    vec output(inputSize);

    it_assert_debug(inputSize % this->size == 0, "Interleaver::interleave(): Invalid length of given vector.");

    for (int i = 0; i < numberOfBlock; i++) {
        for (int j = 0; j < this->size; j++) {
            output[i * this->size + j] = input[i * this->size + this->interleaver[j]];
        }
    }

    return output;
}

bvec Interleaver::deinterleave(const bvec& output)
{
    int outputSize = output.size();
    int numberOfBlock = outputSize / this->size;
    bvec input(outputSize);

    it_assert_debug(outputSize % this->size == 0, "Interleaver::deinterleave(): Invalid length of given vector.");

    for (int i = 0; i < numberOfBlock; i++) {
        for (int j = 0; j < this->size; j++) {
            input[i * this->size + this->interleaver[j]] = output[i * this->size + j];
        }
    }

    return input;
}

vec Interleaver::deinterleave(const vec& output)
{
    int outputSize = output.size();
    int numberOfBlock = outputSize / this->size;
    vec input(outputSize);

    it_assert_debug(outputSize % this->size == 0, "Interleaver::deinterleave(): Invalid length of given vector.");

    for (int i = 0; i < numberOfBlock; i++) {
        for (int j = 0; j < this->size; j++) {
            input[i * this->size + this->interleaver[j]] = output[i * this->size + j];
        }
    }

    return input;
}

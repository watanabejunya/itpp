#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/converters.h>
#include "../../Library/Interleaver/Interleaver.hpp"
#include "../../Library/TurboCode/TurboCode.hpp"
#include "../../Library/TurboCode/PuncturedTurboCode.hpp"

using itpp::bvec;
using itpp::ivec;
using itpp::vec;
using itpp::bin;
using itpp::bmat;
using itpp::mat;
using itpp::round_i;


PuncturedTurboCode::PuncturedTurboCode()
{
    this->puncturePeriod = 0;
}

PuncturedTurboCode::~PuncturedTurboCode()
{
}

bmat PuncturedTurboCode::get_puncture_matrix()
{
    return this->punctureMatrix;
};

int PuncturedTurboCode::get_puncture_period()
{
    return this->puncturePeriod;
};

int PuncturedTurboCode::get_punctured_length()
{
    it_assert_debug(this->puncturePeriod != 0, "PuncturedTurboCode::get_punctured_length(): Punctured_Turbo_Codec: puncture matrix is not set");

    return this->puncturedCodeLength;
};

int PuncturedTurboCode::get_tail_length()
{
    return this->puncturedCodeLength - static_cast<int>(this->informationLength / this->rate);
}

void PuncturedTurboCode::set_puncture_matrix(const bmat& punctureMatrix, bool tailBitUse1, bool tailBitUse2)
{
    int count = 0;
    this->punctureMatrix = punctureMatrix;
    this->puncturePeriod = punctureMatrix.cols();
    this->tailBitUse1 = tailBitUse1;
    this->tailBitUse2 = tailBitUse2;

    it_assert_debug(punctureMatrix.rows() == 3, "PuncturedTurboCode::set_puncture_matrix(): Wrong size of puncture matrix.");
    it_assert_debug(punctureMatrix.cols() != 0, "PuncturedTurboCode::set_puncture_matrix(): Wrong size of puncture matrix.");
    it_assert_debug(this->informationLength % this->puncturePeriod == 0, "PuncturedTurboCode::set_puncture_matrix(): Wrong size of puncture matrix.");

    // Count length of punctured bits.
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < this->puncturePeriod; j++) {
            count += static_cast<int>(this->punctureMatrix(i, j));
        }
    }

    it_assert_debug(count != 0, "PuncturedTurboCode::set_puncture_matrix(): Invalid puncture matrix.");

    this->rate = this->puncturePeriod / (double)(count);

    this->puncturedCodeLength = this->informationLength / this->puncturePeriod * count;
    this->puncturedCodeLength += this->tailBitUse1 ? 2 * (this->constraintLength - 1) : 0;
    this->puncturedCodeLength += this->tailBitUse2 ? 2 * (this->constraintLength - 1) : 0;
}

bvec PuncturedTurboCode::encode(const bvec &input)
{
    bvec output;
    this->encode(input, output);
    return output;
}

void PuncturedTurboCode::encode(const bvec& input, bvec& output)
{
    int index = 0;
    bmat dataMatrix, tailMatrix;
    output.set_size(this->puncturedCodeLength);

    it_assert_debug(this->puncturePeriod != 0, "PuncturedTurboCode::encode(): Puncture matrix is not set.");
    it_assert_debug(input.size() == this->informationLength, "PuncturedTurboCode::encode(): The size of given vector doesn't equal to information length.");

    this->block_encode(input, dataMatrix, tailMatrix);

    // Puncture data bits.
    for (int i = 0; i < this->informationLength; i++) {
        for (int j = 0; j < 3; j++) {
            if (this->punctureMatrix(j, i % this->puncturePeriod) == bin(1)) {
                output[index++] = dataMatrix(j, i);
            }
        }
    }

    // Retrieve the first tail bits.
    if (this->tailBitUse1) {
        for (int i = 0; i < this->constraintLength - 1; i++) {
            output[index++] = tailMatrix(1, i);
            output[index++] = dataMatrix(1, this->informationLength + i);
        }
    }

    // Retrieve the second tail bits.
    if (this->tailBitUse2) {
        for (int i = 0; i < this->constraintLength - 1; i++) {
            output[index++] = tailMatrix(2, i);
            output[index++] = dataMatrix(2, this->informationLength + i);
        }
    }
}

bvec PuncturedTurboCode::decode(const bvec& input)
{
    it_error("PuncturedTurboCode::decode(): Hard-decision decoding is not implemented.");
    return bvec();
}

void PuncturedTurboCode::decode(const bvec& input, bvec& output)
{
    it_error("PuncturedTurboCode::decode(): Hard-decision decoding is not implemented.");
}

bvec PuncturedTurboCode::decode(const vec& input)
{
    bvec output;
    this->decode(input, output);
    return output;
}

void PuncturedTurboCode::decode(const vec& input, bvec& output)
{
    int index = 0;
    mat dataMatrix, tailMatrix;
    dataMatrix.set_size(3, this->informationLength + this->constraintLength - 1);
    tailMatrix.set_size(3, this->constraintLength - 1);

    it_assert_debug(this->puncturePeriod != 0, "PuncturedTurboCode::decode(): Puncture matrix is not set.");
    it_assert_debug(input.size() == this->puncturedCodeLength, "PuncturedTurboCode::decode(): The size of given vector doesn't equal to code length.");

    // Puncture data bits.
    for (int i = 0; i < this->informationLength; i++) {
        for (int j = 0; j < 3; j++) {
            if (this->punctureMatrix(j, i % this->puncturePeriod) == bin(1)) {
                dataMatrix(j, i) = input[index++];
            } else {
                // Insert dummy LLR with same contribution for 0 and 1.
                dataMatrix(j, i) = 0;
            }
        }
    }

    // Retrieve the first tail bits.
    for (int i = 0; i < this->constraintLength - 1; i++) {
        if (this->tailBitUse1) {
            tailMatrix(1, i) = input[index++];
            dataMatrix(1, this->informationLength + i) = input[index++];
        } else {
            // Insert dummy LLR with same contribution for 0 and 1.
            tailMatrix(1, i) = 0;
            dataMatrix(1, this->informationLength + i) = 0;
        }
    }

    // Retrieve the second tail bits.
    for (int i = 0; i < this->constraintLength - 1; i++) {
        if (this->tailBitUse2) {
            tailMatrix(2, i) = input[index++];
            dataMatrix(2, this->informationLength + i) = input[index++];
        } else {
            // Insert dummy LLR with same contribution for 0 and 1.
            tailMatrix(2, i) = 0;
            dataMatrix(2, this->informationLength + i) = 0;
        }
    }

    this->block_decode(dataMatrix, tailMatrix, output);
}

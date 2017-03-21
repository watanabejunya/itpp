#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/comm/rec_syst_conv_code.h>
#include "../../Library/Interleaver/Interleaver.hpp"
#include "../../Library/TurboCode/TurboCode.hpp"
#include "../../Library/Common/helper.hpp"

using std::string;
using itpp::bin;
using itpp::bvec;
using itpp::ivec;
using itpp::vec;
using itpp::bmat;
using itpp::mat;
using itpp::zeros;
using helper::hard_decide_LLRs;


TurboCode::TurboCode()
{
}

TurboCode::~TurboCode()
{
}

void TurboCode::set_parameters(ivec polynomial1, ivec polynomial2, int constraintLength, const ivec& interleaver, int iterations, const string& method)
{
    it_assert_debug(polynomial1.size() == 2, "TurboCode::set_parameters(): This class is supported to 1/3 rate turbo code.");
    it_assert_debug(polynomial2.size() == 2, "TurboCode::set_parameters(): This class is supported to 1/3 rate turbo code.");

    this->encoder1.set_generator_polynomials(polynomial1, constraintLength);
    // this->encoder1.set_method_type(method);
    this->encoder2.set_generator_polynomials(polynomial2, constraintLength);
    // this->encoder2.set_method_type(method);
    this->constraintLength = constraintLength;
    this->rate = 1.0 / 3.0;
    this->set_interleaver(interleaver);
    this->iterations = iterations;
    this->set_decoding_method(method);
}

void TurboCode::set_interleaver(const itpp::ivec& interleaver)
{
    this->informationLength = interleaver.size();
    this->codeLength = informationLength * 3 + (this->constraintLength - 1) * 4;

    this->interleaver.set_interleaver(interleaver);
}

void TurboCode::set_decoding_method(const string& method)
{
    it_assert_debug(method == "LOGMAX" || method == "LOGMAP" || method == "MAP" || method == "TABLE",
        "TurboCode::set_decoding_method(): Invalid decoding method.");

    this->decodingMethod = method;
}

int TurboCode::get_code_length() const
{
    return codeLength;
}

int TurboCode::get_infomation_length() const
{
    return informationLength;
}

double TurboCode::get_rate() const
{
    return this->rate;
}

bvec TurboCode::encode(const itpp::bvec& input)
{
    bvec output;
    this->encode(input, output);
    return output;
}

void TurboCode::encode(const bvec& input, bvec& output)
{
    int index = 0;
    bmat dataMatrix, tailMatrix;
    output.set_size(this->codeLength);

    // Encode and get matrix-format output.
    this->block_encode(input, dataMatrix, tailMatrix);

    // Retrieve the data part.
    for(int i = 0; i < this->informationLength; i++) {
        for (int j = 0; j < 3; j++) {
            output[index++] = dataMatrix(j, i);
        }
    }

    // Retrieve the first tail part.
    for(int i = 0; i < this->constraintLength - 1; i++) {
        output[index++] = tailMatrix(1, i);
        output[index++] = dataMatrix(1, this->informationLength + i);
    }

    // Retrieve the second tail part.
    for(int i = 0; i < this->constraintLength - 1; i++) {
        output[index++] = tailMatrix(2, i);
        output[index++] = dataMatrix(2, this->informationLength + i);
    }
}

bvec TurboCode::decode(const bvec& input)
{
    it_error("TurboCode::decode(): Hard-decision decoding is not implemented.");
    return bvec();
}

void TurboCode::decode(const bvec& input, bvec& output)
{
    it_error("TurboCode::decode(): Hard-decision decoding is not implemented.");
}

bvec TurboCode::decode(const vec& input)
{
    bvec output;
    this->decode(input, output);
    return output;
}

void TurboCode::decode(const vec& input, bvec& output)
{
    int index = 0;
    mat dataMatrix, tailMatrix;
    dataMatrix.set_size(3, this->informationLength + this->constraintLength - 1);
    tailMatrix.set_size(3, this->constraintLength - 1);

    it_assert_debug(input.size() == this->codeLength, "TurboCode::decode(): The size of given vector doesn't equal to code length.");

    // Retrieve the data part.
    for(int i = 0; i < this->informationLength; i++) {
        for (int j = 0; j < 3; j++) {
            dataMatrix(j, i) = input[index++];
        }
    }

    // Retrieve the first tail part.
    for(int i = 0; i < this->constraintLength - 1; i++) {
        tailMatrix(1, i) = input[index++];
        dataMatrix(1, this->informationLength + i) = input[index++];
    }

    // Retrieve the second tail part.
    for(int i = 0; i < this->constraintLength - 1; i++) {
        tailMatrix(2, i) = input[index++];
        dataMatrix(2, this->informationLength + i) = input[index++];
    }

    this->block_decode(dataMatrix, tailMatrix, output);
}

void TurboCode::block_encode(const bvec& input, bmat& dataMatrix, bmat& tailMatrix)
{
    bvec tailBits1, tailBits2;
    bmat dataMatrix1, dataMatrix2;
    dataMatrix.set_size(3, this->informationLength + this->constraintLength - 1);
    tailMatrix.set_size(3, this->constraintLength - 1);

    it_assert_debug(input.length() == this->informationLength, "TurboCode::block_encode(): Parameter error in infomation length.");

    // Encode by each enocder.
    this->encoder1.encode_tail(input, tailBits1, dataMatrix1);
    this->encoder2.encode_tail(this->interleaver.interleave(input), tailBits2, dataMatrix2);

    // Format the data part.
    dataMatrix.set_row(0, input);
    dataMatrix.set_row(1, dataMatrix1.get_col(0));
    dataMatrix.set_row(2, dataMatrix2.get_col(0));

    // Format the tail part.
    tailMatrix.set_row(1, tailBits1);
    tailMatrix.set_row(2, tailBits2);
}

void TurboCode::block_decode(const mat& dataMatrix, const mat& tailMatrix, bvec& output)
{
    vec input, interlaevedInput;
    vec parityLLRs1, parityLLRs2, LLRs1, LLRs2, extrinsicLLRs;
    vec LLRs;

    // Set sizes of vectors.
    input.set_size(this->informationLength + this->constraintLength - 1);
    interlaevedInput.set_size(this->informationLength + this->constraintLength - 1);
    parityLLRs1.set_size(this->informationLength + this->constraintLength - 1);
    parityLLRs2.set_size(this->informationLength + this->constraintLength - 1);
    LLRs1.set_size(this->informationLength);
    LLRs2.set_size(this->informationLength);
    extrinsicLLRs.set_size(this->informationLength);

    // Format the data LLRs.
    input.set_subvector(0, dataMatrix.get_row(0).left(this->informationLength));
    interlaevedInput.set_subvector(0, this->interleaver.interleave(input.left(this->informationLength)));
    parityLLRs1.set_subvector(0, dataMatrix.get_row(1).left(this->informationLength));
    parityLLRs2.set_subvector(0, dataMatrix.get_row(2).left(this->informationLength));

    // Format tail LLRs.
    for(int k = 0; k < this->constraintLength - 1; k++) {
        input[this->informationLength + k] = tailMatrix(1, k);
        interlaevedInput[this->informationLength + k] = tailMatrix(2, k);
        parityLLRs1[this->informationLength + k] = dataMatrix(1, this->informationLength + k);
        parityLLRs2[this->informationLength + k] = dataMatrix(2, this->informationLength + k);
    }

    // Init extrinsic LLRs as 0s.
    extrinsicLLRs.zeros();

    // Start interative decoding.
    for(int j = 0; j < iterations; j++) {
        this->encoder1.log_decode_n2(input, parityLLRs1, extrinsicLLRs, LLRs1, true, this->decodingMethod);
        extrinsicLLRs = this->interleaver.interleave(LLRs1.left(this->informationLength));

        this->encoder2.log_decode_n2(interlaevedInput, parityLLRs2, extrinsicLLRs, LLRs2, true, this->decodingMethod);
        extrinsicLLRs = this->interleaver.deinterleave(LLRs2.left(this->informationLength));
    }

    // Take final bit decisions.
    LLRs = input.left(this->informationLength) + LLRs1.left(this->informationLength) + this->interleaver.deinterleave(LLRs2.left(this->informationLength));
    output = hard_decide_LLRs(LLRs);
}

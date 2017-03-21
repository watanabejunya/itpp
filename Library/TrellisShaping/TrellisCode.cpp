#include <itpp/itbase.h>
#include "../../Library/Common/helper.hpp"
#include "../../Library/TrellisShaping/TrellisCode.hpp"

using itpp::bvec;
using itpp::ivec;
using itpp::vec;
using itpp::cvec;
using itpp::imat;
using itpp::dec2bin;
using itpp::levels2bits;
using itpp::zeros;
using itpp::min_index;
using helper::INFTY;


const static imat generatorStateTable(
    "0 2;"
    "0 2;"
    "1 3;"
    "1 3"
);

const static imat generatorOutputTable(
    "0 3;"
    "3 0;"
    "1 2;"
    "2 1"
);

const static imat inverseSyndromeStateTable(
    "0 1;"
    "0 1"
);

const static imat inverseSyndromeOutputTable(
    "0 1;"
    "3 2"
);

const static imat parityCheckerStateTable(
    "0 2 8 10;"
    "0 2 8 10;"
    "1 3 9 11;"
    "1 3 9 11;"
    "0 2 8 10;"
    "0 2 8 10;"
    "1 3 9 11;"
    "1 3 9 11;"
    "4 6 12 14;"
    "4 6 12 14;"
    "5 7 13 15;"
    "5 7 13 15;"
    "4 6 12 14;"
    "4 6 12 14;"
    "5 7 13 15;"
    "5 7 13 15"
);

const static imat parityCheckerOutputTable(
    "0 1 1 0;"
    "1 0 0 1;"
    "0 1 1 0;"
    "1 0 0 1;"
    "1 0 0 1;"
    "0 1 1 0;"
    "1 0 0 1;"
    "0 1 1 0;"
    "1 0 0 1;"
    "0 1 1 0;"
    "1 0 0 1;"
    "0 1 1 0;"
    "0 1 1 0;"
    "1 0 0 1;"
    "0 1 1 0;"
    "1 0 0 1"
);

const static imat compoundDecoderneStateTable(
    "0 4 6 2;"
    "2 6 4 0;"
    "5 1 3 7;"
    "7 3 1 5;"
    "2 6 4 0;"
    "0 4 6 2;"
    "7 3 1 5;"
    "5 1 3 7"
);

const static imat compoundDecoderOutputTable(
    "0 1 1 0;"
    "0 1 1 0;"
    "1 0 0 1;"
    "1 0 0 1;"
    "0 1 1 0;"
    "0 1 1 0;"
    "1 0 0 1;"
    "1 0 0 1"
);

TrellisCode::TrellisCode()
{
}

TrellisCode::~TrellisCode()
{
}

void TrellisCode::set_type(const TrellisCodeType codeType)
{
    switch (codeType) {
        case TrellisCodeType::Generator:
            this->memorySize = 2;
            this->numberOfStates = 4;
            this->numberOfTransits = 2;
            this->inputBlockLength = 1;
            this->outputBlockLength = 2;
            this->stateTable = generatorStateTable;
            this->outputTable = generatorOutputTable;
            break;
        case TrellisCodeType::InverseSyndrome:
            this->memorySize = 1;
            this->numberOfStates = 2;
            this->numberOfTransits = 2;
            this->inputBlockLength = 1;
            this->outputBlockLength = 2;
            this->stateTable = inverseSyndromeStateTable;
            this->outputTable = inverseSyndromeOutputTable;
            break;
        case TrellisCodeType::ParityChecker:
            this->memorySize = 4;
            this->numberOfStates = 16;
            this->numberOfTransits = 4;
            this->inputBlockLength = 2;
            this->outputBlockLength = 1;
            this->stateTable = parityCheckerStateTable;
            this->outputTable = parityCheckerOutputTable;
            break;
        case TrellisCodeType::CompoundDecoder:
            this->memorySize = 3;
            this->numberOfStates = 8;
            this->numberOfTransits = 4;
            this->inputBlockLength = 2;
            this->outputBlockLength = 1;
            this->stateTable = compoundDecoderneStateTable;
            this->outputTable = compoundDecoderOutputTable;
            break;
        default:
            it_error("TrellisCode::set_type(): Invalid type of code.");
    }

    this->codeType = codeType;
    this->isValid = true;
}

void TrellisCode::init_trellis(const int dataSize)
{
    this->pathMemories.set_size(dataSize + 1, this->numberOfStates);
    this->pathMemories.set_submatrix(0, dataSize, 0, this->numberOfStates - 1, -1);
    this->pathMemories(0, 0) = 0;

    this->outputSymbols.set_size(dataSize + 1, this->numberOfStates);

    this->nodeMetrics.set_size(dataSize + 1, this->numberOfStates);
    this->nodeMetrics.set_submatrix(0, dataSize, 0, this->numberOfStates - 1, INFTY);

    this->forwardMetrics.set_size(dataSize + 1, this->numberOfStates);
    this->forwardMetrics.set_submatrix(0, dataSize, 0, this->numberOfStates - 1, -INFTY);

    this->backwardMetrics.set_size(dataSize + 1, this->numberOfStates);
    this->backwardMetrics.set_submatrix(0, dataSize, 0, this->numberOfStates - 1, -INFTY);

    this->autocorrelation.set_size(dataSize + 1, this->numberOfStates);
    this->autocorrelation.set_submatrix(0, dataSize, 0, this->numberOfStates - 1, to_cvec(zeros(dataSize)));

    this->isReady = true;
}

int TrellisCode::get_memory_size()
{
    return this->memorySize;
}

int TrellisCode::get_number_of_states()
{
    return this->numberOfStates;
}

int TrellisCode::get_number_of_transits()
{
    return this->numberOfTransits;
}

int TrellisCode::get_input_block_length()
{
    return this->inputBlockLength;
}

int TrellisCode::get_output_block_length()
{
    return this->outputBlockLength;
}

ivec TrellisCode::get_previous_states(const int nextState)
{
    it_assert_debug(this->isValid, "TrellisCode::get_previous_states(): Type is not set.");
    it_assert_debug(nextState < this->numberOfStates, "TrellisCode::get_previous_states(): Over the number of state.");

    ivec states;
    for (int currentState = 0; currentState < this->numberOfStates; currentState++) {
        for (int input = 0; input < this->numberOfTransits; input++) {
            if (this->stateTable(currentState, input) == nextState) {
                states.ins(states.size(), currentState);
                break;
            }
        }
    }
    return states;
}

ivec TrellisCode::get_next_states(const int currentState)
{
    it_assert_debug(this->isValid, "TrellisCode::get_next_states(): Type is not set.");
    it_assert_debug(currentState < this->numberOfStates, "TrellisCode::get_next_states(): Over the number of state.");

    return this->stateTable.get_row(currentState);
}

int TrellisCode::get_input(const int currentState, const int nextState)
{
    it_assert_debug(this->isValid, "TrellisCode::get_input(): Type is not set.");
    it_assert_debug(0 <= currentState && currentState < this->numberOfStates, "TrellisCode::get_input(): Over the number of state.");
    it_assert_debug(0 <= nextState && nextState < this->numberOfStates, "TrellisCode::get_input(): Over the number of state.");

    for (int input = 0; input < this->numberOfTransits; input++) {
        if (this->stateTable(currentState, input) == nextState) {
            return input;
        }
    }

    it_error("TrellisCode::get_input(): There is no path between these states.");
    return 0;
}

int TrellisCode::get_next_state(const int currentState, const int input)
{
    it_assert_debug(this->isValid, "TrellisCode::get_next_state(): Type is not set.");
    it_assert_debug(0 <= currentState && currentState < this->numberOfStates, "TrellisCode::get_next_state(): Over the number of currentState.");
    it_assert_debug(0 <= input && input < this->numberOfTransits, "TrellisCode::get_next_state(): Over the number of input.");

    return this->stateTable(currentState, input);
}

int TrellisCode::get_output(const int currentState, const int nextState)
{
    it_assert_debug(this->isValid, "TrellisCode::get_output(): Type is not set.");
    it_assert_debug(0 <= currentState && currentState < this->numberOfStates, "TrellisCode::get_output(): Over the number of state.");
    it_assert_debug(0 <= nextState && nextState < this->numberOfStates, "TrellisCode::get_output(): Over the number of state.");

    for (int input = 0; input < this->numberOfTransits; input++) {
        if (this->stateTable(currentState, input) == nextState) {
            return this->outputTable(currentState, input);
        }
    }

    it_error("TrellisCode::get_output(): There is no path between these states.");
    return 0;
}

bvec TrellisCode::encode(const itpp::bvec& inputs, int state)
{
    int dataSize = inputs.size() / this->inputBlockLength;
    bvec outputs(dataSize * this->outputBlockLength);

    it_assert_debug(this->isValid, "TrellisCode::encode_unterminated(): Type is not set.");
    it_assert_debug(inputs.size() > 0, "TrellisCode::encode_unterminated(): Invalid size of given vector.");
    it_assert_debug(inputs.size() % this->inputBlockLength == 0, "TrellisCode::encode_unterminated(): Invalid size of given vector.");
    it_assert_debug(0 <= state && state < this->numberOfStates, "TrellisCode::encode_unterminated(): Over the number of state.");

    for (int i = 0; i < dataSize; i++) {
        // Retrieve and convert bits to integer value.
        int input = bin2dec(inputs.mid(i * this->inputBlockLength, this->inputBlockLength));

        bvec output = dec2bin(this->outputBlockLength, this->outputTable(state, input));
        outputs.set_subvector(i * this->outputBlockLength, output);

        state = this->stateTable(state, input);
    }
    return outputs;
}

int TrellisCode::trace_back_state(const int node, const int state)
{
    it_assert_debug(this->isValid, "TrellisCode::get_previous_state(): Type is not set.");
    it_assert_debug(0 <= state && state < this->numberOfStates, "TrellisCode::get_previous_state(): Over the number of state.");
    it_assert_debug(this->pathMemories(node, state) >= 0, "TrellisCode::get_previous_state(): Previous path is not set on this node.");

    return this->pathMemories(node, state);
}

cvec TrellisCode::get_symbols(int node, int state)
{
    it_assert_debug(this->isValid, "TrellisCode::get_symbols(): Type is not set.");
    it_assert_debug(0 <= state && state < this->numberOfStates, "TrellisCode::get_symbols(): Over the number of state.");
    it_assert_debug(this->pathMemories(node, state) >= 0, "TrellisCode::get_symbols(): Symbols sequence is not set on this node.");

    cvec symbols(node);
    for (; node > 0; node--) {
        symbols[node - 1] = this->outputSymbols(node, state);
        state = this->trace_back_state(node, state);
    }
    return symbols;
}

int TrellisCode::find_optimum_state(const int node)
{
    return min_index(this->nodeMetrics.get_row(node));
}

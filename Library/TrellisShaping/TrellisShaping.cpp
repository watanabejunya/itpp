#include <itpp/itbase.h>
#include "../../Library/Common/helper.hpp"
#include "../../Library/TrellisShaping/TrellisShaping.hpp"
#include "../../Library/TrellisShaping/TrellisCode.hpp"
#include "../../Library/TrellisShaping/QAM.hpp"
#include "../../Library/ClippingAndFiltering/ClippingAndFiltering.hpp"

using std::max;
using itpp::sqr;
using itpp::pow2;
using itpp::dec2bin;
using itpp::linspace;
using itpp::round_i;
using itpp::bin;
using itpp::vec;
using itpp::bvec;
using itpp::ivec;
using itpp::cvec;
using itpp::mat;
using helper::is_power_of_4;
using helper::demultiplex;
using helper::multiplex;
using helper::hard_decide_LLRs;
using helper::INFTY;
using helper::is_infinity;
using complex = std::complex<double>;


const static mat clippingRatioTable(
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0   0.1   0.2   0.3   0.6;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0   0.5   0.9   1.3   1.4;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0   0.9   1.3   1.7;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0;"
    "-1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0"
);

TrellisShaping::TrellisShaping(const int m, const TSMetricType metric, const TSMappingType mappingType)
{
    this->metric = metric;
    this->rate = (log2(m) - 1.0) / log2(m);

    this->set_constellation(m, mappingType);
    this->set_shaping_method(metric);
    this->set_convolutional_code();

    this->isValid = true;
}

TrellisShaping::~TrellisShaping()
{
}

double TrellisShaping::get_rate()
{
    return this->rate;
}

double TrellisShaping::get_optimum_clipping_ratio(const double scalingFactor)
{
    double clippingRatio = clippingRatioTable(this->QAMModulator.get_k(), round_i(10 * scalingFactor));

    it_assert_debug(clippingRatio > 0, "TrellisShaping::get_optimum_clipping_ratio(): Undefined parameters set is designated.");

    return clippingRatio;
}

void TrellisShaping::set_shaping_method(TSMetricType metricType)
{
    switch (metricType) {
        case TSMetricType::Autocorrelation:
            this->shapingMethod = TSShapingMethod::Sign;
            this->isSetUp = true;
            break;
        case TSMetricType::Clipping:
            this->shapingMethod = TSShapingMethod::Tail;
            break;
        default:
            it_error("Invalid shaping method is given.");
    };
}

void TrellisShaping::set_convolutional_code()
{
    this->generator.set_type(TrellisCodeType::Generator);
    this->inverseSymdrome.set_type(TrellisCodeType::InverseSyndrome);
    this->parityChecker.set_type(TrellisCodeType::ParityChecker);
    this->compoundDecoder.set_type(TrellisCodeType::CompoundDecoder);

    this->numberOfShapingBits = 2;
}

void TrellisShaping::set_constellation(const int m, const TSMappingType mappingType)
{
    it_assert_debug(is_power_of_4(m) && m > 4, "TrellisShaping::set_constellation(): Invalid size of constellation.");

    this->mappingType = mappingType;
    this->QAMModulator.set_constellation(m, mappingType);

    if (this->metric == TSMetricType::Autocorrelation) {
        this->QAMModulator.make_correlation_table();
    }
}


bvec TrellisShaping::add_shaping_bits(const bvec& codeBits, const bvec& shapingBits)
{
    it_assert_debug(codeBits.size() == this->QAMModulator.get_k(), "TrellisShaping::add_shaping_bits(): Invalid code length of given vector.");
    it_assert_debug(shapingBits.size() == this->numberOfShapingBits, "TrellisShaping::add_shaping_bits(): Invalid code length of given vector.");

    bvec shapedBits(codeBits);
    if (this->shapingMethod == TSShapingMethod::Sign) {
        bvec MSBs = codeBits.left(this->numberOfShapingBits) + shapingBits;
        shapedBits.set_subvector(0, MSBs);
    } else if (this->shapingMethod == TSShapingMethod::Tail) {
        bvec LSBs = codeBits.right(this->numberOfShapingBits) + shapingBits;
        shapedBits.set_subvector(shapedBits.size() - LSBs.size(), LSBs);
    }
    return shapedBits;
}

double TrellisShaping::calculate_metric(const cvec& symbols, const cvec& autocorrelation, const int stage)
{
    double branchMetric = 0.0;

    switch (this->metric) {
        case TSMetricType::Autocorrelation:
            for (int m = 1; m < stage; m++) {
                branchMetric += 2.0 * real(conj(autocorrelation[m]) * this->QAMModulator.get_correlation(symbols[stage], symbols[stage - m]));
            }
            if (this->mappingType == TSMappingType::Type2) {
                for (int m = 1; m <= stage; m++) {
                    branchMetric += this->QAMModulator.get_squared_correlation(symbols[stage], symbols[stage - m]);
                }
            }
            return round_i(branchMetric * pow(10.0, 5.0)) * pow(10.0, -5.0);
        case TSMetricType::Clipping:
            branchMetric = norm(symbols[stage] - this->referenceSymbols[stage]);
            return branchMetric;
        default:
            it_error("Invalid metric type is given.");
            return branchMetric;
    }
}


double TrellisShaping::calculate_transit_probability(const int currentState, int nextState,
    const complex& symbol, const double N0, const complex& channel)
{
    const int constellationSize = this->QAMModulator.get_M();
    cvec constellations(constellationSize / 2);
    bvec input = dec2bin(this->numberOfShapingBits, this->compoundDecoder.get_input(currentState, nextState));

    double probability = 1;
    for (int i = 0; i < this->numberOfShapingBits; i++) {
        int index;
        if (this->shapingMethod == TSShapingMethod::Sign) {
            index = i;
        } else {
            index = this->QAMModulator.get_k() - this->numberOfShapingBits + i;
        }

        constellations.set_subvector(0, this->QAMModulator.get_constellation_points(index, input[i]));

        double sum = 0;
        for (int j = 0; j < constellations.size(); j++) {
            sum += (1.0 / constellationSize) / sqrt(M_PI * N0) * exp(- sqr(symbol - channel * constellations[j]) / N0);
        }
        probability *= sum;
    }

    return log(probability);
}

void TrellisShaping::make_reference_symbols(const bvec& codes, double scalingFactor, double clippingRatio)
{
    const int constellationSize = this->QAMModulator.get_M() / 4;
    const int blockLength = this->QAMModulator.get_k();
    clippingRatio = clippingRatio > 0 ? clippingRatio : this->get_optimum_clipping_ratio(scalingFactor);
    ClippingAndFiltering CAFProccesser(clippingRatio);
    bvec MSBs, LSBs;

    it_assert_debug(this->isValid, "TrellisShaping::make_reference_symbols(): Constellation and type are not set.");
    it_assert_debug(this->metric == TSMetricType::Clipping, "TrellisShaping::make_reference_symbols(): This is only for clipping metric.");
    it_assert_debug(codes.size() % blockLength == 0, "TrellisShaping::make_reference_symbols(): Invalid data size of given vecotr.");

    demultiplex(codes, MSBs, LSBs, blockLength - this->numberOfShapingBits, this->numberOfShapingBits);

    // Save original QAM constellation temporarily.
    this->QAMModulator.push_constellation();

    // Set a new QAM constellation for reference symbols.
    scalingFactor *= sqrt(pow2(this->numberOfShapingBits) * (constellationSize - 1) / (constellationSize * 4 - 1));
    this->QAMModulator.set_constellation(constellationSize, this->mappingType, scalingFactor);

    // Make the replica symbols sequence by clipping and set as reference symbols.
    this->referenceSymbols = CAFProccesser.clip_and_filter(this->QAMModulator.modulate_bits(MSBs));
    this->referenceSymbols = CAFProccesser.normalize_attenuation(this->referenceSymbols);

    this->QAMModulator.pop_constellation();
    this->isSetUp = true;
}

bvec TrellisShaping::encode(const bvec& input)
{
    int numberOfNonshapingBits = this->QAMModulator.get_k() - this->numberOfShapingBits;
    bvec MSBs, LSBs;
    bvec output;

    it_assert_debug(input.size() % (this->QAMModulator.get_k() - 1) == 0, "TrellisShaping::encode(): Invalid data size of given vecotr.");

    if (this->shapingMethod == TSShapingMethod::Sign) {
        demultiplex(input, MSBs, LSBs, 1, numberOfNonshapingBits);

        // Encode with inverse symdrome.
        MSBs = this->inverseSymdrome.encode(MSBs);

        multiplex(MSBs, LSBs, output, this->numberOfShapingBits, numberOfNonshapingBits);
    } else {
        demultiplex(input, MSBs, LSBs, numberOfNonshapingBits, 1);

        // Encode with inverse symdrome.
        LSBs = this->inverseSymdrome.encode(LSBs);

        multiplex(MSBs, LSBs, output, numberOfNonshapingBits, this->numberOfShapingBits);
    }

    it_assert_debug(output.size() % this->QAMModulator.get_k() == 0, "TrellisShaping::encode(): Invalid data size of given vecotr.");

    return output;
}

bvec TrellisShaping::decode(const bvec& output)
{
    int numberOfNonshapingBits = this->QAMModulator.get_k() - this->numberOfShapingBits;
    bvec input;
    bvec MSBs, LSBs;

    it_assert_debug(output.size() % this->QAMModulator.get_k() == 0, "TrellisShaping::decode(): Invalid data size of given vecotr.");

    if (this->shapingMethod == TSShapingMethod::Sign) {
        demultiplex(output, MSBs, LSBs, this->numberOfShapingBits, numberOfNonshapingBits);

        // Decode with inverse symdrome.
        MSBs = this->parityChecker.encode(MSBs);

        multiplex(MSBs, LSBs, input, 1, numberOfNonshapingBits);
    } else {
        demultiplex(output, MSBs, LSBs, numberOfNonshapingBits, this->numberOfShapingBits);

        // Decode with inverse symdrome.
        LSBs = this->parityChecker.encode(LSBs);

        multiplex(MSBs, LSBs, input, numberOfNonshapingBits, 1);
    }

    it_assert_debug(input.size() % (this->QAMModulator.get_k() - 1) == 0, "TrellisShaping::decode(): Invalid data size of given vecotr.");

    return input;
}

cvec TrellisShaping::modulate(const bvec& codes)
{
    int blockLength = this->QAMModulator.get_k();
    int dataSize = codes.size() / blockLength;
    bvec tempOriginalBits(blockLength);
    bvec tempShapingBits(this->numberOfShapingBits);
    cvec tempSymbols;

    it_assert_debug(this->isValid, "TrellisShaping::encode_unterminated(): Constellation and type are not set.");
    it_assert_debug(this->isSetUp, "TrellisShaping::encode_unterminated(): Setup has still not been completed.");

    this->generator.init_trellis(dataSize);
    this->generator.nodeMetrics(0, 0) = 0;

    for (int node = 0; node < dataSize; node++) {
        tempSymbols.set_size(node + 1);

        for (int currentState = 0; currentState < this->generator.get_number_of_states(); currentState++) {
            if (is_infinity(this->generator.nodeMetrics(node, currentState))) {
                continue;
            }

            // Get temporarily symbols sequence.
            tempSymbols.set_subvector(0, this->generator.get_symbols(node, currentState));

            ivec nextStates = this->generator.get_next_states(currentState);

            for (int i = 0; i < nextStates.size(); i++) {
                // Make and append a new shped symbol.
                tempOriginalBits.set_subvector(0, codes.mid(node * blockLength, blockLength));
                tempShapingBits.set_subvector(0, dec2bin(this->numberOfShapingBits, this->generator.get_output(currentState, nextStates[i])));
                tempSymbols[node] = this->QAMModulator.modulate_bits(this->add_shaping_bits(tempOriginalBits, tempShapingBits)).get(0);

                // Calculate branch metric.
                double branchMetric = this->calculate_metric(tempSymbols, this->generator.autocorrelation(node, currentState), node);

                if (this->generator.nodeMetrics(node, currentState) + branchMetric < this->generator.nodeMetrics(node + 1, nextStates[i])) {
                    // Update the metric.
                    this->generator.nodeMetrics(node + 1, nextStates[i]) = round(pow(10.0, 5.0) * (this->generator.nodeMetrics(node, currentState) + branchMetric)) * pow(10.0, -5.0);

                    // Reserve previous state and selected outputs.
                    this->generator.pathMemories(node + 1, nextStates[i]) = currentState;
                    this->generator.outputSymbols(node + 1, nextStates[i]) = tempSymbols[node];
                }
            }
        }

        if (this->metric == TSMetricType::Autocorrelation) {
            // Update autocorrelation for next stage.
            for (int currentState = 0; currentState < this->generator.get_number_of_states(); currentState++) {
                if (is_infinity(this->generator.nodeMetrics(node + 1, currentState))) {
                    continue;
                }

                int previousState = this->generator.pathMemories(node + 1, currentState);

                // Get temporarily symbols sequence.
                tempSymbols.set_subvector(0, this->generator.get_symbols(node + 1, currentState));

                for (int m = 1; m <= node; m++) {
                    complex correlation = this->generator.autocorrelation(node, previousState)[m]
                        + this->QAMModulator.get_correlation(tempSymbols[node] , tempSymbols[node - m]);

                    double real = round_i(correlation.real() * pow(10.0, 5.0)) * pow(10.0, -5.0);
                    double imag = round_i(correlation.imag() * pow(10.0, 5.0)) * pow(10.0, -5.0);

                    this->generator.autocorrelation(node + 1, currentState)[m] = complex(real, imag);
                }
            }
        }

    }

    // Find minimum metric.
    int state = this->generator.find_optimum_state(dataSize);

    // Trace back to calculate the output symbol.
    return this->generator.get_symbols(dataSize, state);
}

bvec TrellisShaping::demodulate_bits(const cvec& symbols, const double N0)
{
    bvec codes;

    it_assert_debug(this->isValid, "TrellisShaping::demodulate_bits(): Constellation and type are not set.");
    it_assert_debug(this->isSetUp, "TrellisShaping::demodulate_bits(): Setup has still not been completed.");

    codes = hard_decide_LLRs(this->QAMModulator.demodulate_soft_bits(symbols, N0));

    return this->decode(codes);
}

bvec TrellisShaping::demodulate_bits(const cvec& symbols, const cvec& channel, const double N0)
{
    bvec codes;

    it_assert_debug(this->isValid, "TrellisShaping::demodulate_bits(): Constellation and type are not set.");
    it_assert_debug(this->isSetUp, "TrellisShaping::demodulate_bits(): Setup has still not been completed.");

    codes = hard_decide_LLRs(this->QAMModulator.demodulate_soft_bits(symbols, channel, N0));

    return this->decode(codes);
}

vec TrellisShaping::demodulate_soft_bits(const cvec& symbols, const double N0)
{
    int dataSize = symbols.size();
    int blockLength = this->QAMModulator.get_k();
    int numberOfNonshapingBits = blockLength - this->numberOfShapingBits;
    vec LLRs, nonshapedLLRs, shapedLLRs;

    it_assert_debug(this->isValid, "TrellisShaping::demodulate_soft_bits(): Constellation and type are not set.");
    it_assert_debug(this->isSetUp, "TrellisShaping::demodulate_soft_bits(): Setup has still not been completed.");

    // Caluclulate raw LLRs for non-shaped bits.
    LLRs = this->QAMModulator.demodulate_soft_bits(symbols, N0);

    if (this->shapingMethod == TSShapingMethod::Sign) {
        demultiplex(LLRs, shapedLLRs, nonshapedLLRs, this->numberOfShapingBits, numberOfNonshapingBits);
    } else {
        demultiplex(LLRs, nonshapedLLRs, shapedLLRs, numberOfNonshapingBits, this->numberOfShapingBits);
    }

    // Start max-log-MAP algorithm for shaped bits.
    this->compoundDecoder.init_trellis(dataSize);
    this->compoundDecoder.forwardMetrics(0, 0) = log(1);

    // First, calculate joint probabilities for forward recursion.
    for (int node = 1; node < dataSize; node++) {
        for (int currentState = 0; currentState < this->compoundDecoder.get_number_of_states(); currentState++) {
            ivec previousStates = this->compoundDecoder.get_previous_states(currentState);

            for (int i = 0; i < previousStates.size(); i++) {
                if (is_infinity(this->compoundDecoder.forwardMetrics(node-1, previousStates[i]))) {
                    continue;
                }

                double metric = this->calculate_transit_probability(previousStates[i], currentState, symbols[node-1], N0)
                    + this->compoundDecoder.forwardMetrics(node-1, previousStates[i]);

                this->compoundDecoder.forwardMetrics(node, currentState) = max(this->compoundDecoder.forwardMetrics(node, currentState), metric);
            }
        }
    }

    // Next, obtain LLR of each shaped bit.
    shapedLLRs.set_size(dataSize * 1);

    for (int node = 0; node < dataSize; node++) {
        double metricFor0 = - INFTY;
        double metricFor1 = - INFTY;

        for (int currentState = 0; currentState < this->compoundDecoder.get_number_of_states(); currentState++) {
            ivec nextStates = this->compoundDecoder.get_next_states(currentState);

            for (int i = 0; i < nextStates.size(); i++) {
                int output = this->compoundDecoder.get_output(currentState, nextStates[i]);

                double metric = this->calculate_transit_probability(currentState, nextStates[i], symbols[node], N0)
                    + this->compoundDecoder.forwardMetrics(node, currentState);

                if (output == 0) {
                    metricFor0 = max(metricFor0, metric);
                } else if (output == 1) {
                    metricFor1 = max(metricFor1, metric);
                }
            }
        }
        shapedLLRs[node] = metricFor0 - metricFor1;
    }

    if (this->shapingMethod == TSShapingMethod::Sign) {
        multiplex(shapedLLRs, nonshapedLLRs, LLRs, 1, numberOfNonshapingBits);
    } else {
        multiplex(nonshapedLLRs, shapedLLRs, LLRs, numberOfNonshapingBits, 1);
    }

    return LLRs;
}

vec TrellisShaping::demodulate_soft_bits(const cvec& symbols, const cvec& channel, const double N0)
{
    int dataSize = symbols.size();
    int blockLength = this->QAMModulator.get_k();
    int numberOfNonshapingBits = blockLength - this->numberOfShapingBits;
    vec LLRs, nonshapedLLRs, shapedLLRs;

    it_assert_debug(this->isValid, "TrellisShaping::demodulate_soft_bits(): Constellation and type are not set.");
    it_assert_debug(this->isSetUp, "TrellisShaping::demodulate_soft_bits(): Setup has still not been completed.");
    it_assert_debug(channel.size() == dataSize, "TrellisShaping::demodulate_soft_bits(): Invalid size of given vector.");

    // Caluclulate raw LLRs for non-shaped bits.
    LLRs = this->QAMModulator.demodulate_soft_bits(symbols, channel, N0);

    if (this->shapingMethod == TSShapingMethod::Sign) {
        demultiplex(LLRs, shapedLLRs, nonshapedLLRs, this->numberOfShapingBits, numberOfNonshapingBits);
    } else {
        demultiplex(LLRs, nonshapedLLRs, shapedLLRs, numberOfNonshapingBits, this->numberOfShapingBits);
    }

    // Start max-log-MAP algorithm for shaped bits.
    this->compoundDecoder.init_trellis(dataSize);
    this->compoundDecoder.forwardMetrics(0, 0) = log(1);

    // First, calculate joint probabilities for forward recursion.
    for (int node = 1; node < dataSize; node++) {
        for (int currentState = 0; currentState < this->compoundDecoder.get_number_of_states(); currentState++) {
            ivec previousStates = this->compoundDecoder.get_previous_states(currentState);

            for (int i = 0; i < previousStates.size(); i++) {
                if (is_infinity(this->compoundDecoder.forwardMetrics(node-1, previousStates[i]))) {
                    continue;
                }

                double metric = this->calculate_transit_probability(previousStates[i], currentState, symbols[node-1], N0, channel[node-1])
                    + this->compoundDecoder.forwardMetrics(node-1, previousStates[i]);

                this->compoundDecoder.forwardMetrics(node, currentState)
                    = max(this->compoundDecoder.forwardMetrics(node, currentState), metric);
            }
        }
    }

    // Next, obtain LLR of each shaped bit.
    shapedLLRs.set_size(dataSize * 1);

    for (int node = 0; node < dataSize; node++) {
        double metricFor0 = - INFTY;
        double metricFor1 = - INFTY;

        for (int currentState = 0; currentState < this->compoundDecoder.get_number_of_states(); currentState++) {
            ivec nextStates = this->compoundDecoder.get_next_states(currentState);

            for (int i = 0; i < nextStates.size(); i++) {
                int output = this->compoundDecoder.get_output(currentState, nextStates[i]);

                double metric = this->calculate_transit_probability(currentState, nextStates[i], symbols[node], N0, channel[node])
                    + this->compoundDecoder.forwardMetrics(node, currentState);

                if (output == 0) {
                    metricFor0 = max(metricFor0, metric);
                } else if (output == 1) {
                    metricFor1 = max(metricFor1, metric);
                }
            }
        }
        shapedLLRs[node] = metricFor0 - metricFor1;
    }

    if (this->shapingMethod == TSShapingMethod::Sign) {
        multiplex(shapedLLRs, nonshapedLLRs, LLRs, 1, numberOfNonshapingBits);
    } else {
        multiplex(nonshapedLLRs, shapedLLRs, LLRs, numberOfNonshapingBits, 1);
    }

    return LLRs;
}

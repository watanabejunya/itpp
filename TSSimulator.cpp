#include <iostream>
#include <iomanip>
#include <itpp/itbase.h>
#include "config.hpp"
#include "./Library/Common/helper.hpp"
#include "./Library/Counter/Counters.hpp"
#include "./Library/Channel/Channels.hpp"
#include "./Library/Interleaver/Interleavers.hpp"
#include "./Library/Histogram/Histogram.hpp"
#include "./Library/OFDM/OFDM.hpp"
#include "./Library/TrellisShaping/TrellisShaping.hpp"
#include "./Library/TurboCode/PuncturedTurboCode.hpp"


using std::string;
using svec = std::vector<std::string>;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::ofstream;
using std::ios;
using itpp::GlobalRNG_randomize;
using itpp::vec;
using itpp::bvec;
using itpp::ivec;
using itpp::cvec;
using itpp::bmat;
using itpp::cmat;
using itpp::randb;
using itpp::linspace_fixed_step;
using itpp::pow2i;
using itpp::floor_i;
using itpp::dB;
using itpp::inv_dB;
using itpp::Qfunc;
using helper::parse_arguments;
using helper::join;
using helper::assume_PAPR;
using helper::get_average_power;
using helper::hard_decide_LLRs;


void calculate_instantaneous_CCDF()
{
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    Histogram PSDHistogram(0, 13, 1001);
    OFDM OFDM(NUMBER_OF_SUBCARRIERS, 0, OVERSAMPLING_FACTOR);
    TrellisShaping TS(CONSTELLATION_SIZE, TSMetricType::Autocorrelation, TSMappingType(TS_MAPPING_TYPE));

    ofstream file(join("./Result/TS/CCDF_Norm_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE , ".dat"), ios::out);

    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cout << "trial: " << i + 1 << "\r";

        transmittedBits = randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);

        transmittedCodes = TS.encode(transmittedBits);

        transmittedSymbols = TS.modulate(transmittedCodes);

        transmittedSignals = OFDM.modulate(transmittedSymbols);

        PSDHistogram.read_normalized_power(transmittedSignals);
    }

    cout << "\n";
    PSDHistogram.print_CCDF(file);
}

void calculate_average_power_reduction_capability()
{
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    double capability = 0;
    QAM QAM(CONSTELLATION_SIZE);
    TrellisShaping TS(CONSTELLATION_SIZE, TSMetricType::Autocorrelation, TSMappingType(TS_MAPPING_TYPE));
    OFDM OFDM(NUMBER_OF_SUBCARRIERS, GUARD_INTERVAL, OVERSAMPLING_FACTOR);

    ofstream file(join("./Result/TS/APR_Capa_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE , ".dat"), ios::out);

    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cout << "trial: " << i + 1 << "\r";

        transmittedBits = randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);

        transmittedCodes = TS.encode(transmittedBits);

        transmittedSymbols = TS.modulate(transmittedCodes);

        transmittedSignals = OFDM.modulate(transmittedSymbols);

        capability += dB(1.0 / get_average_power(transmittedSignals)) / (double)(NUMBER_OF_TRIALS);
    }

    cout << capability << "[dB]" << "\n";
    file << capability << endl;
}

void calculate_uncoded_BER_over_AWGN()
{
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    cvec receivedSignals;
    cvec receivedSymbols;
    vec receivedLLRs;
    bvec receivedBits;
    TrellisShaping TS(CONSTELLATION_SIZE, TSMetricType::Autocorrelation, TSMappingType(TS_MAPPING_TYPE));
    OFDM OFDM(NUMBER_OF_SUBCARRIERS);
    AWGNChannel AWGNCannel;
    BERCounter BERCounter;
    SERCounter SERCounter{QAM(CONSTELLATION_SIZE)};
    vec EbN0s = linspace_fixed_step(-5.0, 25.0);
    vec SNRs = TS.get_rate() * log2(CONSTELLATION_SIZE) * inv_dB(EbN0s);
    ofstream file1(join("./Result/TS/BER_Uncod_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE , ".dat"), ios::out);
    ofstream file2(join("./Result/TS/SER_Uncod_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE , ".dat"), ios::out);

    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();
        SERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);

            transmittedCodes = TS.encode(transmittedBits);
            transmittedSymbols = TS.modulate(transmittedCodes);

            transmittedSignals = OFDM.modulate(transmittedSymbols);

            AWGNCannel.set_noise(get_average_power(transmittedSymbols) / SNRs[i]);
            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = OFDM.demodulate(receivedSignals);

            receivedBits = TS.demodulate_bits(receivedSymbols, AWGNCannel.get_noise());

            BERCounter.count(transmittedBits, receivedBits);
            SERCounter.count(transmittedSymbols, receivedSymbols, AWGNCannel.get_noise());
        }

        cout << "\n";
        file1 << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;
        file2 << dB(SNRs[i]) << " " << SERCounter.get_errorrate() << endl;

        if (BERCounter.get_errorrate() < BER_LIMIT) break;
    }
}

void calculate_uncoded_BER_over_Rayleigh()
{
    bvec transmittedBits, receivedBits;
    bvec transmittedCodes;
    vec receivedLLRs;
    cvec transmittedSymbols, receivedSymbols;
    cvec transmittedSignals, receivedSignals;
    OFDM OFDM(NUMBER_OF_SUBCARRIERS, GUARD_INTERVAL);
    TrellisShaping TS(CONSTELLATION_SIZE, TSMetricType::Autocorrelation, TSMappingType(TS_MAPPING_TYPE));
    AWGNChannel AWGNCannel;
    StaticFadingChannel fadingChannel(NUMBER_OF_TAPS, DelayProfileType::Uniform);
    fadingChannel.set_fft_size(NUMBER_OF_SUBCARRIERS);
    double averagePower;
    BERCounter BERCounter;
    SERCounter SERCounter{QAM(CONSTELLATION_SIZE)};
    vec EbN0s = linspace_fixed_step(5.0, 25.0);
    vec SNRs = TS.get_rate() * log2(CONSTELLATION_SIZE) * inv_dB(EbN0s);
    ofstream file1(join("./Result/TS/BER_Uncod_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);
    ofstream file2(join("./Result/TS/SER_Uncod_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);


    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();
        SERCounter.clear();

        averagePower = 0.0;
        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            transmittedBits = randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);

            transmittedCodes = TS.encode(transmittedBits);
            transmittedSymbols = TS.modulate(transmittedCodes);

            averagePower += get_average_power(transmittedSymbols) / NUMBER_OF_TRIALS;
        }
        AWGNCannel.set_noise(averagePower / SNRs[i]);

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);

            transmittedCodes = TS.encode(transmittedBits);
            transmittedSymbols = TS.modulate(transmittedCodes);

            transmittedSignals = OFDM.modulate(transmittedSymbols);

            fadingChannel.init();
            receivedSignals = fadingChannel(transmittedSignals);

            receivedSignals = AWGNCannel(receivedSignals);

            receivedSymbols = OFDM.demodulate(receivedSignals);

            receivedBits = TS.demodulate_bits(receivedSymbols, fadingChannel.get_channel_cofficients(), AWGNCannel.get_noise());

            BERCounter.count(transmittedBits, receivedBits);
            SERCounter.count(transmittedSymbols, receivedSymbols, fadingChannel.get_channel_cofficients(), AWGNCannel.get_noise());
        }

        cout << "\n";
        file1 << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;
        file2 << dB(SNRs[i]) << " " << SERCounter.get_errorrate() << endl;

        if (BERCounter.get_errorrate() < BER_LIMIT) break;
    }
}

void calculate_turbo_coded_BER_over_AWGN()
{
    const int tailSize = 12;
    const int informationLength = (NUMBER_OF_BITS - 2) * 124;
    const int codeLength = informationLength * (NUMBER_OF_BITS - 1) / (NUMBER_OF_BITS - 2) + tailSize;
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    cvec receivedSignals;
    cvec receivedSymbols;
    vec receivedLLRs;
    bvec receivedBits;
    TrellisShaping TS(CONSTELLATION_SIZE, TSMetricType::Autocorrelation, TSMappingType(TS_MAPPING_TYPE));
    OFDM OFDM(NUMBER_OF_SUBCARRIERS);
    RandomInterleaver randomInterleaver((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);
    ivec generator("13 15");
    bmat punctureMatrix(3, 2 * (NUMBER_OF_BITS - 2));
    punctureMatrix.zeros();
    punctureMatrix.set_submatrix(0, 0, 0, 2 * (NUMBER_OF_BITS - 2) - 1, 1);
    punctureMatrix(1, 1) = 1;
    punctureMatrix(2, 1) = 1;
    PuncturedTurboCode turboCode;
    turboCode.set_parameters(generator, generator, 4, SRandomInterleaver::make_interleaver(informationLength), NUMBER_OF_TURBO_ITERATIONS, "LOGMAP");
    turboCode.set_puncture_matrix(punctureMatrix);
    AWGNChannel AWGNCannel;
    BERCounter BERCounter;
    FERCounter FERCounter(informationLength);
    vec EbN0s = linspace_fixed_step(6.0, 23.0, 0.5);
    vec SNRs = TS.get_rate() * turboCode.get_rate() * log2(CONSTELLATION_SIZE) * inv_dB(EbN0s);
    ofstream file1(join("./Result/TS/BER_Turbo_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);
    ofstream file2(join("./Result/TS/FER_Turbo_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);

    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();
        FERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(informationLength);

            // Encode by turbo code and fill random bits up to N * (B - 1).
            transmittedCodes = turboCode.encode(transmittedBits);
            transmittedCodes.ins(transmittedCodes.size(), randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS - codeLength));

            // Random interleave the bit sequence to avoid burst errors.
            transmittedCodes = randomInterleaver.interleave(transmittedCodes);

            // Proccess trellis shaping to reduce PAPR.
            transmittedCodes = TS.encode(transmittedCodes);
            transmittedSymbols = TS.modulate(transmittedCodes);

            transmittedSignals = OFDM.modulate(transmittedSymbols);

            // Transmit signals over AWGN channel.
            AWGNCannel.set_noise(get_average_power(transmittedSymbols) / SNRs[i]);
            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = OFDM.demodulate(receivedSignals);

            // Demodulate and calculate LLRs with trellis.
            receivedLLRs = TS.demodulate_soft_bits(receivedSymbols, AWGNCannel.get_noise());

            // Random interleave the bit sequence to avoid burst errors.
            receivedLLRs = randomInterleaver.deinterleave(receivedLLRs);

            // Decode information bits by turbo decoder.
            receivedBits = turboCode.decode(receivedLLRs.left(codeLength));

            BERCounter.count(transmittedBits, receivedBits);
            FERCounter.count(transmittedBits, receivedBits);
        }

        cout << "\n";
        file1 << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;
        file2 << dB(SNRs[i]) << " " << FERCounter.get_errorrate() << endl;

        if (BERCounter.get_errorrate() < BER_LIMIT || FERCounter.get_errorrate() < FER_LIMIT) break;
    }
}

void calculate_turbo_coded_BER_over_Rayleigh()
{
    const int tailSize = 12;
    const int informationLength = (NUMBER_OF_BITS - 2) * 124;
    const int codeLength = informationLength * (NUMBER_OF_BITS - 1) / (NUMBER_OF_BITS - 2) + tailSize;
    bvec transmittedBits, receivedBits;
    bvec transmittedCodes;
    vec receivedLLRs;
    cvec transmittedSymbols, receivedSymbols;
    cvec transmittedSignals, receivedSignals;
    OFDM OFDM(NUMBER_OF_SUBCARRIERS, GUARD_INTERVAL);
    TrellisShaping TS(CONSTELLATION_SIZE, TSMetricType::Autocorrelation, TSMappingType(TS_MAPPING_TYPE));
    RandomInterleaver randomInterleaver((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);
    ivec generator("13 15");
    bmat punctureMatrix(3, 2 * (NUMBER_OF_BITS - 2));
    PuncturedTurboCode turboCode;
    punctureMatrix.zeros();
    punctureMatrix.set_submatrix(0, 0, 0, 2 * (NUMBER_OF_BITS - 2) - 1, 1);
    punctureMatrix(1, 1) = 1;
    punctureMatrix(2, 6) = 1;
    turboCode.set_parameters(generator, generator, 4, SRandomInterleaver::make_interleaver(informationLength), NUMBER_OF_TURBO_ITERATIONS, "LOGMAP");
    turboCode.set_puncture_matrix(punctureMatrix);
    AWGNChannel AWGNCannel;
    StaticFadingChannel fadingChannel(NUMBER_OF_TAPS, DelayProfileType::Uniform);
    fadingChannel.set_fft_size(NUMBER_OF_SUBCARRIERS);
    BERCounter BERCounter;
    FERCounter FERCounter(informationLength);
    double averagePower;
    vec EbN0s = linspace_fixed_step(10.0, 25.0);
    vec SNRs = TS.get_rate() * turboCode.get_rate() * NUMBER_OF_BITS * inv_dB(EbN0s);
    ofstream file1(join("./Result/TS/BER_Turbo_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);
    ofstream file2(join("./Result/TS/FER_Turbo_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);


    averagePower = 0.0;
    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cerr << "trial: " << setw(5) << i + 1 << " Average Power (dB): " << setw(2) << averagePower * NUMBER_OF_TRIALS / i << "\r";

        transmittedBits = randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);

        transmittedCodes = TS.encode(transmittedBits);
        transmittedSymbols = TS.modulate(transmittedCodes);

        averagePower += get_average_power(transmittedSymbols) / NUMBER_OF_TRIALS;
    }

    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        AWGNCannel.set_noise(averagePower / SNRs[i]);
        BERCounter.clear();
        FERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " SNR[dB]: " << setw(2) << dB(SNRs[i]) << " FER: " << FERCounter.get_errorrate() << "\r";

            transmittedBits = randb(informationLength);

            // Encode by turbo code and fill random bits up to N * (B - 1).
            transmittedCodes = turboCode.encode(transmittedBits);
            transmittedCodes.ins(transmittedCodes.size(), randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS - codeLength));

            // Random interleave the bit sequence to avoid burst errors.
            transmittedCodes = randomInterleaver.interleave(transmittedCodes);

            transmittedCodes = TS.encode(transmittedCodes);
            transmittedSymbols = TS.modulate(transmittedCodes);

            transmittedSignals = OFDM.modulate(transmittedSymbols);

            fadingChannel.init();
            receivedSignals = fadingChannel(transmittedSignals);

            receivedSignals = AWGNCannel(receivedSignals);

            receivedSymbols = OFDM.demodulate(receivedSignals);

            // Demodulate and calculate LLRs with trellis.
            receivedLLRs = TS.demodulate_soft_bits(receivedSymbols, fadingChannel.get_channel_cofficients(), AWGNCannel.get_noise());

            // Random interleave the bit sequence to avoid burst errors.
            receivedLLRs = randomInterleaver.deinterleave(receivedLLRs);

            // Decode information bits by turbo decoder.
            receivedBits = turboCode.decode(receivedLLRs.left(codeLength));

            BERCounter.count(transmittedBits, receivedBits);
            FERCounter.count(transmittedBits, receivedBits);
        }

        cout << "\n";
        file1 << EbN0s[i] << " " << BERCounter.get_errorrate() << endl;
        file2 << dB(SNRs[i]) << " " << FERCounter.get_errorrate() << endl;

        if (BERCounter.get_errorrate() < BER_LIMIT || FERCounter.get_errorrate() < FER_LIMIT) break;
    }
}

void calculate_theoretical_BER_over_AWGN()
{
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    QAM QAM(CONSTELLATION_SIZE);
    TrellisShaping TS(CONSTELLATION_SIZE, TSMetricType::Autocorrelation, TSMappingType(TS_MAPPING_TYPE));
    Histogram symbolHistogram(0, sqrt(CONSTELLATION_SIZE) - 1, sqrt(CONSTELLATION_SIZE));
    vec SNRs = inv_dB(linspace_fixed_step(0.0, 30.0));
    double capability = 0.0;
    vec symbolDistribution;
    double occurrenceProbability;
    double numberOfMSBError, numberOfLSBError, numberOfBitError;
    double BER;
    double SER;
    ofstream file1(join("./Result/TS/BER_Theory_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);
    ofstream file2(join("./Result/TS/SER_Theory_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);

    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cerr << "trial: " << setw(5) << i + 1 << "\r";

        transmittedBits = randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);

        transmittedCodes = TS.encode(transmittedBits);
        transmittedSymbols = TS.modulate(transmittedCodes);

        capability += 1.0 / get_average_power(transmittedSymbols) / NUMBER_OF_TRIALS;

        symbolHistogram.read_number(QAM.demodulate_real(transmittedSymbols));
    }

    // Calculate the symbol occurrence probability.
    symbolDistribution = symbolHistogram.get_discrete_histogram();

    occurrenceProbability = symbolDistribution[0];
    for (int i = 1; i < symbolDistribution.size() - 1; i++) {
        occurrenceProbability += 2.0 * symbolDistribution[i];
    }
    occurrenceProbability += symbolDistribution[symbolDistribution.size() - 1];

    numberOfMSBError = 1.0 / 2.0;
    numberOfLSBError = 1.0 / 2.0;
    numberOfBitError = numberOfMSBError + numberOfLSBError * 2.5;

    for (int i = 0; i < SNRs.size(); i++) {
        SER = occurrenceProbability * Qfunc(sqrt(3.0 * capability / (CONSTELLATION_SIZE - 1.0) * SNRs[i])) * (2.0 - occurrenceProbability * Qfunc(sqrt(3.0 * capability / (CONSTELLATION_SIZE - 1.0) * SNRs[i])));
        BER = numberOfBitError / (NUMBER_OF_BITS - 1) * SER;

        file1 << dB(SNRs[i]) << " " << BER << endl;
        file2 << dB(SNRs[i]) << " " << SER << endl;

        cout << " Es/N0(dB): " << setw(2) << dB(SNRs[i]) << " SER: " << SER << setw(2) << " BER: " << BER << "\n";

        if (BER < BER_LIMIT || SER < SER_LIMIT) break;
    }
}

void calculate_theoretical_BER_over_Rayleigh()
{
    bvec transmittedBits;
    bvec transmittedCodes;
    cvec transmittedSymbols;
    QAM QAM(CONSTELLATION_SIZE);
    TrellisShaping TS(CONSTELLATION_SIZE, TSMetricType::Autocorrelation, TSMappingType(TS_MAPPING_TYPE));
    Histogram symbolHistogram(0, sqrt(CONSTELLATION_SIZE) - 1, sqrt(CONSTELLATION_SIZE));
    vec SNRs = inv_dB(linspace_fixed_step(10.0, 50.0));
    double capability = 0.0;
    vec symbolDistribution;
    double occurrenceProbability;
    double numberOfMSBError, numberOfLSBError, numberOfBitError;
    double BER;
    double SER;
    ofstream file1(join("./Result/TS/BER_Theory_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);
    ofstream file2(join("./Result/TS/SER_Theory_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", "_Type-", TS_MAPPING_TYPE, ".dat"), ios::out);

    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cerr << "trial: " << setw(5) << i + 1 << "\r";

        transmittedBits = randb((NUMBER_OF_BITS - 1) * NUMBER_OF_SUBCARRIERS);

        transmittedCodes = TS.encode(transmittedBits);
        transmittedSymbols = TS.modulate(transmittedCodes);

        capability += 1.0 / get_average_power(transmittedSymbols) / NUMBER_OF_TRIALS;

        symbolHistogram.read_number(QAM.demodulate_real(transmittedSymbols));
    }

    // Calculate the symbol occurrence probability.
    symbolDistribution = symbolHistogram.get_discrete_histogram();

    occurrenceProbability = symbolDistribution[0];
    for (int i = 1; i < symbolDistribution.size() - 1; i++) {
        occurrenceProbability += 2.0 * symbolDistribution[i];
    }
    occurrenceProbability += symbolDistribution[symbolDistribution.size() - 1];

    // Calucalate the number of bit errors.
    numberOfMSBError = 1.0 / 2.0;
    numberOfLSBError = 1.0 / 2.0;
    numberOfBitError = numberOfMSBError + numberOfLSBError * 2.5;

    for (int i = 0, size = SNRs.size(); i < size; i++) {
        SER = occurrenceProbability * (1.0 - pow(2.0 * (CONSTELLATION_SIZE - 1.0) / (3.0 * capability * SNRs[i]) + 1.0, -1.0))
            - pow(occurrenceProbability, 2.0) / 4.0 * (pow(108.0 * capability * SNRs[i] / (CONSTELLATION_SIZE - 1.0) + 36.0, -1.0)
            + pow(16.0 * capability * SNRs[i] / (CONSTELLATION_SIZE - 1.0) + 4.0, -1.0)
            + pow(21.0 * capability * SNRs[i] / (CONSTELLATION_SIZE - 1.0) + 6.0, -1.0));
        BER = numberOfBitError / (NUMBER_OF_BITS - 1) * SER;

        file1 << dB(SNRs[i]) << " " << BER << endl;
        file2 << dB(SNRs[i]) << " " << SER << endl;

        cout << " Es/N0(dB): " << setw(2) << dB(SNRs[i]) << " SER: " << SER << setw(2) << " BER: " << BER << "\n";

        if (BER < BER_LIMIT || SER < SER_LIMIT) break;
    }
}


int main(const int argc, const char* argv[])
{
    svec arguments = parse_arguments(argc, argv);
    GlobalRNG_randomize();

    if (arguments.size() <= 1) {
        cerr << "No arguments error." << endl;
        exit(1);
    }

    if (arguments[1] == "normal") {
        calculate_instantaneous_CCDF();
    } else if (arguments[1] == "average") {
        calculate_average_power_reduction_capability();
    } else if (arguments[1] == "uncoded") {
        if (arguments[2] == "awgn") {
            calculate_uncoded_BER_over_AWGN();
        } else if (arguments[2] == "fading") {
            calculate_uncoded_BER_over_Rayleigh();
        } else {
            cerr << "Invalid arguments error." << endl;
            exit(1);
        }
    } else if (arguments[1] == "turbo") {
        if (arguments[2] == "awgn") {
            calculate_turbo_coded_BER_over_AWGN();
        } else if (arguments[2] == "fading") {
            calculate_turbo_coded_BER_over_Rayleigh();
        } else {
            cerr << "Invalid arguments error." << endl;
            exit(1);
        }
    } else if (arguments[1] == "theory") {
        if (arguments[2] == "awgn") {
            calculate_theoretical_BER_over_AWGN();
        } else if (arguments[2] == "fading") {
            calculate_theoretical_BER_over_Rayleigh();
        } else {
            cerr << "Invalid arguments error." << endl;
            exit(1);
        }
    } else {
        cerr << "Invalid arguments error." << endl;
        exit(1);
    }

    return 0;
}

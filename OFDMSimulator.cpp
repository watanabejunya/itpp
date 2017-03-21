#include <iostream>
#include <iomanip>
#include <itpp/comm/modulator.h>
#include "config.hpp"
#include "./Library/Common/helper.hpp"
#include "./Library/Counter/Counters.hpp"
#include "./Library/Channel/Channels.hpp"
#include "./Library/Interleaver/Interleavers.hpp"
#include "./Library/Histogram/Histogram.hpp"
#include "./Library/OFDM/OFDM.hpp"
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
using itpp::ivec;
using itpp::bvec;
using itpp::cvec;
using itpp::bmat;
using itpp::cmat;
using itpp::randb;
using itpp::linspace_fixed_step;
using itpp::pow2;
using itpp::fact;
using itpp::floor_i;
using itpp::dB;
using itpp::inv_dB;
using itpp::Qfunc;
using itpp::to_bvec;
using itpp::QAM;
using helper::parse_arguments;
using helper::join;
using helper::get_average_power;
using helper::hard_decide_LLRs;


void calculate_instantaneous_CCDF()
{
    bvec transmittedBits;
    cvec transmittedSymbols;
    cvec transmittedSignals;
    Histogram PSDHistogram(0, 13, 1001);
    QAM QAM(CONSTELLATION_SIZE);
    OFDM OFDM(NUMBER_OF_SUBCARRIERS, 0, OVERSAMPLING_FACTOR);
    ofstream file(join("./Result/OFDM/CCDF_Norm_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);


    for (int i = 0; i < NUMBER_OF_TRIALS; i++) {
        cout << "trial: " << i + 1 << "\r";

        transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

        transmittedSymbols = QAM.modulate_bits(transmittedBits);

        transmittedSignals = OFDM.modulate(transmittedSymbols);

        PSDHistogram.read_normalized_power(transmittedSignals);
    }

    cout << "\n";
    PSDHistogram.print_CCDF(file);
}

void calculate_uncoded_BER_over_AWGN()
{
    bvec transmittedBits, receivedBits;
    cvec transmittedSymbols, receivedSymbols;
    cvec transmittedSignals, receivedSignals;
    QAM QAM(CONSTELLATION_SIZE);
    OFDM OFDM(NUMBER_OF_SUBCARRIERS);
    AWGNChannel AWGNCannel;
    BERCounter BERCounter;
    SERCounter SERCounter{QAM};
    vec EbN0s = linspace_fixed_step(-3.0, 22.0);
    vec SNRs = NUMBER_OF_BITS * inv_dB(EbN0s);
    ofstream file1(join("./Result/OFDM/BER_Uncod_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);
    ofstream file2(join("./Result/OFDM/SER_Uncod_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);


    for (size_t i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();
        SERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

            transmittedSymbols = QAM.modulate_bits(transmittedBits);

            transmittedSignals = OFDM.modulate(transmittedSymbols);

            AWGNCannel.set_noise(get_average_power(transmittedSignals) / SNRs[i]);
            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = OFDM.demodulate(receivedSignals);

            receivedBits = hard_decide_LLRs(QAM.demodulate_soft_bits(receivedSymbols, AWGNCannel.get_noise()));

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
    cvec transmittedSymbols, receivedSymbols;
    cvec transmittedSignals, receivedSignals;
    QAM QAM(CONSTELLATION_SIZE);
    OFDM OFDM(NUMBER_OF_SUBCARRIERS, GUARD_INTERVAL);
    AWGNChannel AWGNCannel;
    StaticFadingChannel fadingChannel(NUMBER_OF_TAPS, DelayProfileType::Uniform);
    fadingChannel.set_fft_size(NUMBER_OF_SUBCARRIERS);
    BERCounter BERCounter;
    SERCounter SERCounter{QAM};
    vec EbN0s = linspace_fixed_step(0.0, 22.0);
    vec SNRs = NUMBER_OF_BITS * inv_dB(EbN0s);
    ofstream file1(join("./Result/OFDM/BER_Uncod_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);
    ofstream file2(join("./Result/OFDM/SER_Uncod_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);


    for (size_t i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();
        SERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);

            transmittedSymbols = QAM.modulate_bits(transmittedBits);

            transmittedSignals = OFDM.modulate(transmittedSymbols);

            fadingChannel.init();
            receivedSignals = fadingChannel(transmittedSignals);

            AWGNCannel.set_noise(get_average_power(receivedSignals) / SNRs[i]);
            receivedSignals = AWGNCannel(receivedSignals);

            receivedSymbols = OFDM.demodulate(receivedSignals);

            receivedBits = hard_decide_LLRs(QAM.demodulate_soft_bits(receivedSymbols, fadingChannel.get_channel_cofficients(), AWGNCannel.get_noise()));

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
    const int codeLength = informationLength * NUMBER_OF_BITS / (NUMBER_OF_BITS - 2) + tailSize;
    bvec transmittedBits, receivedBits;
    bvec transmittedCodes;
    vec receivedLLRs;
    cvec transmittedSymbols, receivedSymbols;
    cvec transmittedSignals, receivedSignals;
    QAM QAM(CONSTELLATION_SIZE);
    OFDM OFDM(NUMBER_OF_SUBCARRIERS);
    RandomInterleaver interleaver(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);
    ivec generator("13 15");
    bmat punctureMatrix(3, (NUMBER_OF_BITS - 2));
    PuncturedTurboCode turboCode;
    punctureMatrix.zeros();
    punctureMatrix.set_submatrix(0, 0, 0, (NUMBER_OF_BITS - 2) - 1, 1);
    punctureMatrix(1, 1) = 1;
    punctureMatrix(2, (NUMBER_OF_BITS - 2) * 3 / 2 - 6) = 1;
    turboCode.set_parameters(generator, generator, 4, SRandomInterleaver::make_interleaver(informationLength), NUMBER_OF_TURBO_ITERATIONS, "LOGMAP");
    turboCode.set_puncture_matrix(punctureMatrix);
    AWGNChannel AWGNCannel;
    BERCounter BERCounter;
    FERCounter FERCounter(informationLength);
    vec EbN0s = linspace_fixed_step(5.0, 25.0, 0.5);
    vec SNRs = turboCode.get_rate() * NUMBER_OF_BITS * inv_dB(EbN0s);
    ofstream file1(join("./Result/OFDM/BER_Turbo_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);
    ofstream file2(join("./Result/OFDM/FER_Turbo_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);


    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        BERCounter.clear();
        FERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " Eb/N0(dB): " << setw(2) << EbN0s[i] << " BER: " << BERCounter.get_errorrate() << "\r";

            transmittedBits = randb(informationLength);

            // Encode by turbo code and fill random bits up to N * (B - 1).
            transmittedCodes = turboCode.encode(transmittedBits);
            transmittedCodes.ins(transmittedCodes.size(), randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS - codeLength));

            // Random interleave the bit sequence to avoid burst errors.
            transmittedCodes = interleaver.interleave(transmittedCodes);

            transmittedSymbols = QAM.modulate_bits(transmittedCodes);

            transmittedSignals = OFDM.modulate(transmittedSymbols);

            AWGNCannel.set_noise(get_average_power(transmittedSignals) / SNRs[i]);
            receivedSignals = AWGNCannel(transmittedSignals);

            receivedSymbols = OFDM.demodulate(receivedSignals);

            receivedLLRs = QAM.demodulate_soft_bits(receivedSymbols, AWGNCannel.get_noise());

            // Deinterleave the bit sequence to avoid burst errors.
            receivedLLRs = interleaver.deinterleave(receivedLLRs);

            // Demodulate and calculate LLRs with trellis.
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
    const int informationLength = (NUMBER_OF_BITS - 2) * floor_i(NUMBER_OF_SUBCARRIERS - tailSize / (double)(NUMBER_OF_BITS));
    const int codeLength = informationLength * NUMBER_OF_BITS / (NUMBER_OF_BITS - 2) + tailSize;
    bvec transmittedBits, receivedBits;
    bvec transmittedCodes;
    vec receivedLLRs;
    cvec transmittedSymbols, receivedSymbols;
    cvec transmittedSignals, receivedSignals;
    QAM QAM(CONSTELLATION_SIZE);
    OFDM OFDM(NUMBER_OF_SUBCARRIERS, GUARD_INTERVAL);
    RandomInterleaver interleaver(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS);
    ivec generator("13 15");
    bmat punctureMatrix(3, (NUMBER_OF_BITS - 2));
    PuncturedTurboCode turboCode;
    punctureMatrix.zeros();
    punctureMatrix.set_submatrix(0, 0, 0, NUMBER_OF_BITS - 3, 1);
    punctureMatrix(1, 1) = 1;
    punctureMatrix(2, (NUMBER_OF_BITS - 2) * 3 / 2 - 6) = 1;
    turboCode.set_parameters(generator, generator, 4, SRandomInterleaver::make_interleaver(informationLength), NUMBER_OF_TURBO_ITERATIONS, "LOGMAP");
    turboCode.set_puncture_matrix(punctureMatrix);
    AWGNChannel AWGNCannel;
    StaticFadingChannel fadingChannel(NUMBER_OF_TAPS, DelayProfileType::Uniform);
    fadingChannel.set_fft_size(NUMBER_OF_SUBCARRIERS);
    BERCounter BERCounter;
    FERCounter FERCounter(informationLength);
    vec EbN0s = linspace_fixed_step(10.0, 30.0);
    vec SNRs = turboCode.get_rate() * NUMBER_OF_BITS * inv_dB(EbN0s);
    ofstream file1(join("./Result/OFDM/BER_Turbo_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);
    ofstream file2(join("./Result/OFDM/FER_Turbo_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);

    for (int i = 0, size = EbN0s.size(); i < size; i++) {
        AWGNCannel.set_noise(1.0 / SNRs[i]);
        BERCounter.clear();
        FERCounter.clear();

        for (int j = 0; j < NUMBER_OF_TRIALS; j++) {
            cerr << "trial: " << setw(5) << j + 1 << " SNR[dB]: " << setw(2) << dB(SNRs[i]) << " FER: " << FERCounter.get_errorrate() << "\r";

            transmittedBits = randb(informationLength);

            // Encode by turbo code and fill random bits up to N * (B - 1).
            transmittedCodes = turboCode.encode(transmittedBits);
            transmittedCodes.ins(transmittedCodes.size(), randb(NUMBER_OF_BITS * NUMBER_OF_SUBCARRIERS - codeLength));

            // Random interleave the bit sequence to avoid burst errors.
            transmittedCodes = interleaver.interleave(transmittedCodes);

            transmittedSymbols = QAM.modulate_bits(transmittedCodes);

            transmittedSignals = OFDM.modulate(transmittedSymbols);

            fadingChannel.init();
            receivedSignals = fadingChannel(transmittedSignals);

            receivedSignals = AWGNCannel(receivedSignals);

            receivedSymbols = OFDM.demodulate(receivedSignals);

            receivedLLRs = QAM.demodulate_soft_bits(receivedSymbols, fadingChannel.get_channel_cofficients(), AWGNCannel.get_noise());

            // Deinterleave the bit sequence to avoid burst errors.
            receivedLLRs = interleaver.deinterleave(receivedLLRs);

            // Demodulate and calculate LLRs with trellis.
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
    vec SNRs = inv_dB(linspace_fixed_step(0.0, 30.0));
    double occurrenceProbability;
    double BER;
    double SER;
    ofstream file1(join("./Result/OFDM/BER_Theory_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);
    ofstream file2(join("./Result/OFDM/SER_Theory_AWGN_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);

    occurrenceProbability = 2.0 * (sqrt(CONSTELLATION_SIZE) - 1.0) / sqrt(CONSTELLATION_SIZE);

    for (int i = 0; i < SNRs.size(); i++) {
        SER = occurrenceProbability * Qfunc(sqrt(3.0 / (CONSTELLATION_SIZE - 1.0) * SNRs[i]))
            * (2.0 - occurrenceProbability * Qfunc(sqrt(3.0 / (CONSTELLATION_SIZE - 1.0) * SNRs[i])));
        BER = SER * (1.0 / NUMBER_OF_BITS);

        file1 << dB(SNRs[i]) << " " << BER << endl;
        file2 << dB(SNRs[i]) << " " << SER << endl;

        cout << " Es/N0(dB): " << setw(2) << dB(SNRs[i]) << " SER: " << SER << setw(2) << " BER: " << BER << "\n";

        if (BER < BER_LIMIT || SER < SER_LIMIT) break;
    }
}

void calculate_theoretical_BER_over_Rayleigh()
{
    vec SNRs = inv_dB(linspace_fixed_step(0.0, 30.0));
    double occurrenceProbability;
    double BER;
    double SER;
    ofstream file1(join("./Result/OFDM/BER_Theory_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);
    ofstream file2(join("./Result/OFDM/SER_Theory_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);

    occurrenceProbability = 2.0 * (sqrt(CONSTELLATION_SIZE) - 1.0) / sqrt(CONSTELLATION_SIZE);

    for (int i = 0, size = SNRs.size(); i < size; i++) {
        SER = occurrenceProbability * (1.0 - pow(2.0 * (CONSTELLATION_SIZE - 1.0) / (3.0 * SNRs[i]) + 1.0, -1.0))
            - pow(occurrenceProbability, 2.0) / 4.0 * (pow(108.0 * SNRs[i] / (CONSTELLATION_SIZE - 1.0) + 36.0, -1.0)
            + pow(16.0 * SNRs[i] / (CONSTELLATION_SIZE - 1.0) + 4.0, -1.0)
            + pow(21.0 * SNRs[i] / (CONSTELLATION_SIZE - 1.0) + 6.0, -1.0));
        BER = SER * (1.0 / NUMBER_OF_BITS);

        file1 << dB(SNRs[i]) << " " << BER << endl;
        file2 << dB(SNRs[i]) << " " << SER << endl;

        cout << " Es/N0(dB): " << setw(2) << dB(SNRs[i]) << " SER: " << SER << setw(2) << " BER: " << BER << "\n";

        if (BER < BER_LIMIT || SER < SER_LIMIT) break;
    }
}

void calculate_outage_probability_over_Rayleigh()
{
    double dataRate = NUMBER_OF_BITS - 2;
    vec SNRs = inv_dB(linspace_fixed_step(10.0, 25.0));
    ofstream file(join("./Result/OFDM/FER_Outage_Rayleigh_", CONSTELLATION_SIZE, "-QAM_", NUMBER_OF_SUBCARRIERS, "-Sub", ".dat"), ios::out);

    for (int i = 0; i < SNRs.size(); i++) {
        double sum = 0.0;
        for (int l = 0; l < NUMBER_OF_TAPS; l++) {
            sum += pow(NUMBER_OF_TAPS * (pow2(dataRate) - 1.0) / SNRs[i], l) / fact(l);
        }

        double outageProbability = 1.0 - exp(- NUMBER_OF_TAPS * (pow2(dataRate) - 1.0) / SNRs[i]) * sum;

        file << dB(SNRs[i]) << " " << outageProbability << endl;

        cout << " Es/N0(dB): " << setw(2) << dB(SNRs[i]) << " FER: " << outageProbability << "\n";
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
    } else if (arguments[1] == "outage") {
        calculate_outage_probability_over_Rayleigh();
    } else {
        cerr << "Invalid arguments error." << endl;
        exit(1);
    }

    cout << "Simulation finished successfully!" << endl;
    return 0;
}

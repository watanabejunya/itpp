#include <iostream>
#include <complex>
#include <string>
#include <sstream>
#include <itpp/itbase.h>
#include <itpp/signal/transforms.h>
#include "../../Library/Common/helper.hpp"

using std::endl;
using std::ofstream;
using std::abs;
using std::norm;
using std::pow;
using std::sqrt;
using std::string;
using svec = std::vector<std::string>;
using std::stringstream;
using std::move;
using itpp::bin;
using itpp::vec;
using itpp::bvec;
using itpp::cvec;
using itpp::levels2bits;
using itpp::dB;
using itpp::fft;
using itpp::ifft;


namespace helper {
    const double INFTY = 1e30;

    void __join(std::stringstream&)
    {
    }

    bool contains(double value, vec vector)
    {
        for (int i = 0; i < vector.size(); i++) {
            if ((vector[i] - value) < pow(10.0, -5.0)) {
                return true;
            }
        }
        return false;
    }

    svec parse_arguments(const int argc, const char* argv[])
    {
        svec vector;
        for (int i = 0; i < argc; i++) {
            vector.push_back(string(argv[i]));
        }
        return vector;
    }

    bool is_infinity(double number)
    {
        return abs(number) >= INFTY - 1;
    }

    bool is_power_of_2(const int number)
    {
        return (number > 1) && ((number & (number - 1)) == 0);
    }

    bool is_power_of_4(const int number)
    {
        return (number > 1) && ((number & (number - 1)) == 0) && ((levels2bits(number) & 1) == 0);
    }

    int mod(const int number, const int law)
    {
        int reminder = number % law;
        return reminder >= 0 ? reminder : reminder + law;
    }

    void demultiplex(const bvec& inputs, bvec& outputs1, bvec& outputs2, int number1, int number2)
    {
        int blockLength = number1 + number2;
        int dataSize = inputs.size() / blockLength;

        it_assert_debug(inputs.size() % blockLength == 0, "helper::demultiplex(): Data size of gievn vectors dpesn't match.");

        outputs1.set_size(number1 * dataSize);
        outputs2.set_size(number2 * dataSize);
        for (int i = 0; i < dataSize; i++) {
            outputs1.set_subvector(i * number1, inputs.mid(i * blockLength, number1));
            outputs2.set_subvector(i * number2, inputs.mid(i * blockLength + number1, number2));
        }
    }

    void demultiplex(const vec& inputs, vec& outputs1, vec& outputs2, int number1, int number2)
    {
        int blockLength = number1 + number2;
        int dataSize = inputs.size() / blockLength;

        it_assert_debug(inputs.size() % blockLength == 0, "helper::demultiplex(): Data size of gievn vectors dpesn't match.");

        outputs1.set_size(number1 * dataSize);
        outputs2.set_size(number2 * dataSize);
        for (int i = 0; i < dataSize; i++) {
            outputs1.set_subvector(i * number1, inputs.mid(i * blockLength, number1));
            outputs2.set_subvector(i * number2, inputs.mid(i * blockLength + number1, number2));
        }
    }

    void multiplex(const bvec& inputs1, const bvec& inputs2, bvec& outputs, int number1, int number2)
    {
        int blockLength = number1 + number2;
        int dataSize = inputs1.size() / number1;

        it_assert_debug(inputs1.size() / number1 == inputs2.size() / number2, "helper::multiplex(): Data size of gievn vectors dpesn't match.");

        outputs.set_size(blockLength * dataSize);
        for (int i = 0; i < dataSize; i++) {
            outputs.set_subvector(i * blockLength, inputs1.mid(i * number1, number1));
            outputs.set_subvector(i * blockLength + number1, inputs2.mid(i * number2, number2));
        }
    }

    void multiplex(const vec& inputs1, const vec& inputs2, vec& outputs, int number1, int number2)
    {
        int blockLength = number1 + number2;
        int dataSize = inputs1.size() / number1;

        it_assert_debug(inputs1.size() / number1 == inputs2.size() / number2, "helper::multiplex(): Data size of gievn vectors dpesn't match.");

        outputs.set_size(blockLength * dataSize);
        for (int i = 0; i < dataSize; i++) {
            outputs.set_subvector(i * blockLength, inputs1.mid(i * number1, number1));
            outputs.set_subvector(i * blockLength + number1, inputs2.mid(i * number2, number2));
        }
    }

    cvec rotate_right(cvec vec, const int number)
    {
        it_assert_debug(0 <= number && number < vec.size(), "helper::rotate_right(): Invalid number of shifts.");

        vec.shift_right(vec.right(number));
        return vec;
    }

    cvec rotate_left(cvec vec, const int number)
    {
        it_assert_debug(0 <= number && number < vec.size(), "helper::rotate_left(): Invalid number of shifts.");

        vec.shift_left(vec.left(number));
        return vec;
    }

    cvec circular_convolute(const cvec& vec1, const cvec& vec2, int size)
    {
        size = (size == 0) ? vec1.size() : size;
        const int law = vec1.size();
        cvec vec3(size);

        it_assert_debug(vec1.size() == vec2.size(), "helper::circular_convolute(): Different data sizes of two vectors are given.");

        vec3.zeros();
        for (size_t i = 0; i < size; i++) {
            for (size_t j = 0; j < law; j++) {
                vec3[i] += vec1[mod(j, law)] * vec2[mod(i - j, law)];
            }
        }
        return vec3;
    }

    cvec oversample(const cvec& signals1, const int j)
    {
        const int size = signals1.size();
        const int half = size / 2;

        it_assert_debug(is_power_of_2(size), "helper::oversample(): Data size must be power of 2.");
        it_assert_debug(is_power_of_2(j), "helper::oversample(): Data size must be power of 2.");

        return concat(signals1.right(half), itpp::zeros_c(size * (j - 1)), signals1.left(half));
    }

    cvec downsample(const cvec& signals1, const int j)
    {
        const int size = signals1.size() / j;
        const int half = size / 2;

        it_assert_debug(is_power_of_2(size), "helper::oversample(): Data size must be power of 2.");
        it_assert_debug(is_power_of_2(j), "helper::oversample(): Data size must be power of 2.");

        return concat(signals1.right(half), signals1.left(half));
    }

    cvec add_cyclic_prefix(const cvec& signals1, int length)
    {
        it_assert_debug(0 < signals1.size(), "helper::add_cyclic_prefix(): Data size must be positive number.");
        it_assert_debug(0 < length, "helper::add_cyclic_prefix(): Guard interval must be positive number.");

        cvec signals2(signals1.size() + length);
        signals2.set_subvector(0, signals1.right(length));
        signals2.set_subvector(length, signals1);
        return signals2;
    }

    cvec remove_cyclic_prefix(const cvec& signals1, int length)
    {
        it_assert_debug(0 < signals1.size(), "helper::remove_cyclic_prefix(): Data size must be positive number.");
        it_assert_debug(0 < length, "helper::remove_cyclic_prefix(): Guard interval must be positive number.");

        return signals1.right(signals1.size() - length);
    }

    cvec fft_normalized(const cvec& signals)
    {
        return fft(signals) / sqrt(signals.size());
    }

    cvec ifft_normalized(const cvec& signals)
    {
        return ifft(signals) * sqrt(signals.size());
    }

    bvec hard_decide_LLRs(const vec& LLRs)
    {
        int size = LLRs.size();

        it_assert_debug(size, "helper::hard_decide(): Data size too short.");

        bvec bits(size);
        for (int i = 0; i < size; i++) {
            bits[i] = (LLRs[i] > 0) ? bin(0) : bin(1);
        }
        return bits;
    }

    double assume_PAPR(const cvec& signal)
    {
        double peak = 0, average = 0;

        for (size_t i = 0, size = signal.size(); i < size ; i++) {
            average += norm(signal[i]) / (double)size;

            if (norm(signal[i]) > peak) {
                peak = norm(signal[i]);
            }
        }

        // Convert the power ratio X to (10 * log10(X)) dB.
        return dB(peak / average);
    }

    double get_average_power(const itpp::cvec& signal)
    {
        double average = 0;
        for (size_t i = 0, size = signal.size(); i < size ; i++) {
            average += norm(signal[i]) / (double)size;
        }
        return average;
    }
}

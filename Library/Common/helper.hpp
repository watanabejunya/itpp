#ifndef helper_hpp
#define helper_hpp

#include <string>
#include <limits>
#include <sstream>
#include <fstream>
#include <itpp/itbase.h>

namespace helper
{
    extern const double INFTY;

    void __join(std::stringstream&);                                                // Join variables and make string.
    template<typename First, typename... Rest>
    void __join(std::stringstream& stream, const First& head, const Rest&... rest)  // Join variables and make string.
    {
        stream << head;
        __join(stream, rest...);
    }
    template<typename First, typename... Rest>
    std::string join(const First& first, const Rest&... rest)                       // Join variables and make string.
    {
        std::stringstream stream;
        __join(stream, first, rest...);
        return stream.str();
    }
    template<typename First>
    void __dump(const First& head)                                                  // Print put variables and then exit.
    {
        std::cerr << head << std::endl;
    }
    template<typename First, typename... Rest>
    void __dump(const First& head, const Rest&... rest)                             // Print put variables and then exit.
    {
        std::cerr << head << std::endl;
        __dump(rest...);
    }
    template<typename First, typename... Rest>
    std::string dump(const First& first, const Rest&... rest)                       // Print put variables and then exit.
    {
        __dump(first, rest...);
        exit(1);
    }
    template<typename T>
    bool contains(T value, itpp::Vec<T> vector)                                     // Check if the vector contains the value.
    {
        for (int i = 0; i < vector.size(); i++) {
            if (vector[i] == value) {
                return true;
            }
        }
        return false;
    }
    std::vector<std::string> parse_arguments(const int, const char**);              // Parse arguments of simulation as strings.
    bool is_infinity(double);                                                       // Check whether the diven number is maximum.
    bool is_power_of_2(int num);                                                    // Check whether the given number is power of 2.
    bool is_power_of_4(int num);                                                    // Check whether the given number is power of 4.
    int mod(int number, int law);                                                   // Get the positive division reminder.
    void demultiplex(const itpp::bvec&, itpp::bvec&, itpp::bvec&, int, int);        // Multiplex two vectors into one vector.
    void demultiplex(const itpp::vec&, itpp::vec&, itpp::vec&, int, int);           // Multiplex two vectors into one vector.
    void multiplex(const itpp::bvec&, const itpp::bvec&, itpp::bvec&, int, int);    // Demultiplex one vector into two vectors.
    void multiplex(const itpp::vec&, const itpp::vec&, itpp::vec&, int, int);       // Demultiplex one vector into two vectors.
    itpp::cvec rotate_right(itpp::cvec, int = 1);                                   // Circular shift each element to right.
    itpp::cvec rotate_left(itpp::cvec, int = 1);                                    // Circular shift each element to left.
    itpp::cvec circular_convolute(const itpp::cvec&, const itpp::cvec&, int = 0);   // Circular convolute two vectors.
    itpp::cvec oversample(const itpp::cvec&, int);                                  // Oversample the sequence.
    itpp::cvec downsample(const itpp::cvec&, int);                                  // Downsample the sequence.
    itpp::cvec add_cyclic_prefix(const itpp::cvec&, int);                           // Add cyclic prefix to the sequence.
    itpp::cvec remove_cyclic_prefix(const itpp::cvec&, int);                        // Remove cyclic prefix from the sequence.
    itpp::cvec fft_normalized(const itpp::cvec&);                                   // FFT without changing average energy.
    itpp::cvec ifft_normalized(const itpp::cvec&);                                  // IFFT without changing average energy.
    itpp::bvec hard_decide_LLRs(const itpp::vec&);                                  // Hard decide LLR and convert to binary code.
    double assume_PAPR(const itpp::cvec&);                                          // Assume PAPR of the sequene.
    double get_average_power(const itpp::cvec&);                                    // Calculate average power of complex signals.
};

#endif

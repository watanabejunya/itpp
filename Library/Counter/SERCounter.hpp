#ifndef SERCounter_hpp
#define SERCounter_hpp

#include <itpp/itbase.h>
#include <itpp/comm/modulator.h>
#include "../../Library/Counter/ERCounter.hpp"

class SERCounter: public ERCounter
{
    private:
        int blockSize;                                                                  // Number of block sizes.
        itpp::Modulator<std::complex<double>> modulator;                                // Symbol modulator.

    public:
        SERCounter(const itpp::Modulator<std::complex<double>>&);                       // Constructor.
        virtual ~SERCounter();                                                          // Deconstructor.
        void set_modulator(const itpp::Modulator<std::complex<double>>&);               // Set a modulator.
        void count(const itpp::cvec&, const itpp::cvec&, double);                       // Count symbol errors.
        void count(const itpp::cvec&, const itpp::cvec&, const itpp::cvec&, double);    // Count symbol errors.
};

#endif

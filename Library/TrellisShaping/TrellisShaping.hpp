#ifndef TrellisShaping_hpp
#define TrellisShaping_hpp

#include <itpp/itbase.h>
#include "../../Library/TrellisShaping/TrellisCode.hpp"
#include "../../Library/TrellisShaping/QAM.hpp"

enum class TSMetricType {Autocorrelation, Clipping};
enum class TSShapingMethod {Sign, Tail};

class TrellisShaping
{
    private:
        TSMetricType metric;
        TSMappingType mappingType;
        TSShapingMethod shapingMethod;
        int numberOfShapingBits;
        double rate;
        bool isValid = false;
        bool isSetUp = false;
        TrellisCode generator;
        TrellisCode inverseSymdrome;
        TrellisCode parityChecker;
        TrellisCode compoundDecoder;
        QAM QAMModulator;
        itpp::cvec referenceSymbols;

        double get_optimum_clipping_ratio(double);
        void set_constellation(const int, const TSMappingType mappingType);
        void set_shaping_method(const TSMetricType);
        void set_convolutional_code();
        itpp::bvec add_shaping_bits(const itpp::bvec&, const itpp::bvec&);
        double calculate_metric(const itpp::cvec&, const itpp::cvec&, int);
        double calculate_transit_probability(const int, const int,
            const std::complex<double>&, const double, const std::complex<double>& = 1.0);
        itpp::bvec decode(const itpp::bvec&);

    public:
        TrellisShaping(const int m, const TSMetricType, const TSMappingType);
        ~TrellisShaping();
        double get_rate();
        void make_reference_symbols(const itpp::bvec&, double, double = -1);
        itpp::bvec encode(const itpp::bvec&);
        itpp::cvec modulate(const itpp::bvec&);
        itpp::bvec demodulate_bits(const itpp::cvec&, const double);
        itpp::bvec demodulate_bits(const itpp::cvec&, const itpp::cvec&, const double);
        itpp::vec demodulate_soft_bits(const itpp::cvec&, const double);
        itpp::vec demodulate_soft_bits(const itpp::cvec&, const itpp::cvec&, const double);
};

#endif

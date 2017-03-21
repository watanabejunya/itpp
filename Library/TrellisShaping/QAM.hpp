#ifndef QAM_hpp
#define QAM_hpp

#include <itpp/itbase.h>
#include <itpp/comm/modulator.h>

enum class TSMappingType {Type1 = 1, Type2 = 2};

class QAM: public itpp::QAM
{
    private:
        itpp::cvec tempSymbols;
        itpp::ivec tempMappers;
        itpp::cmat correlationTable;
        itpp::mat squaredCorrelationTable;

    public:
        QAM();
        QAM(int);
        ~QAM();
        void set_constellation(const int, const TSMappingType, const double = 1.0);
        void push_constellation();
        void pop_constellation();
        itpp::cvec get_constellation_points(const int index, const itpp::bin& input);
        itpp::cvec get_adjacent_symbols(const int index);
        void make_correlation_table();
        std::complex<double> get_correlation(const std::complex<double>&, const std::complex<double>&);
        double get_squared_correlation(const std::complex<double>&, const std::complex<double>&);
        itpp::ivec demodulate_real(const itpp::cvec&);
};

#endif

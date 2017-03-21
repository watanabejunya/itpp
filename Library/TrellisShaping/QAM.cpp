#include <itpp/itbase.h>
#include <itpp/comm/modulator.h>
#include "../../Library/Common/helper.hpp"
#include "../../Library/TrellisShaping/QAM.hpp"

using itpp::bin;
using itpp::ivec;
using itpp::cvec;
using itpp::linspace;
using itpp::round_i;
using itpp::floor_i;
using itpp::to_ivec;
using itpp::to_cvec;
using helper::is_power_of_4;
using complex = std::complex<double>;

QAM::QAM(): itpp::QAM()
{
}

QAM::QAM(int M): itpp::QAM(M)
{
}

QAM::~QAM()
{
}

void QAM::set_constellation(const int m, const TSMappingType type, const double scalingFactor)
{
    double averagePower = sqrt((m - 1) * 2.0 / 3.0);
    cvec baseSymbols(4);
    baseSymbols[0] = complex(1, 1);
    baseSymbols[1] = complex(-1, 1);
    baseSymbols[2] = complex(1, -1);
    baseSymbols[3] = complex(-1, -1);
    cvec newSymbols = baseSymbols;

    it_assert_debug(is_power_of_4(m), "TS::set_constellation(): Invalid size of constellation.");

    if (type == TSMappingType::Type1) {
        for (int n = 16; n <= m; n *= 4) {
            newSymbols.set_size(n);
            newSymbols.set_subvector(0, baseSymbols + complex(sqrt(n) / 2, sqrt(n) / 2));
            for (int i = 0; i < n / 4; i++) {
                newSymbols[n / 4 + i] = conj(-newSymbols[i]);
                newSymbols[n * 2 / 4 + i] = conj(newSymbols[i]);
                newSymbols[n * 3 / 4 + i] = -newSymbols[i];
            }
            baseSymbols = newSymbols;
        }
    } else {
        for (int n = 16; n <= m; n *= 4) {
            newSymbols.set_size(n);
            for (int i = 0; i < n / 4; i++) {
                newSymbols[i] = baseSymbols[i] + complex(sqrt(n) / 2, sqrt(n) / 2);
                newSymbols[n / 4 + i] = baseSymbols[i] + complex(-sqrt(n) / 2, sqrt(n) / 2);
                newSymbols[n * 2 / 4 + i] = baseSymbols[i] + complex(sqrt(n) / 2, -sqrt(n) / 2);
                newSymbols[n * 3 / 4 + i] = baseSymbols[i] + complex(-sqrt(n) / 2, -sqrt(n) / 2);
            }
            baseSymbols = newSymbols;
        }
    }

    newSymbols *= scalingFactor / averagePower;

    this->set(newSymbols, to_ivec(linspace(0, m - 1, m)));
}

void QAM::push_constellation()
{
    this->tempSymbols = this->get_symbols();
    this->tempMappers = this->get_bits2symbols();
}

void QAM::pop_constellation()
{
    this->set(this->tempSymbols, this->tempMappers);
    this->tempSymbols.set_size(0);
    this->tempMappers.set_size(0);
}

cvec QAM::get_constellation_points(const int index, const bin& input)
{
    if (input == bin(0)) {
         return this->modulate(this->S0.get_row(index));
    }
    return this->modulate(this->S1.get_row(index));
}

void QAM::make_correlation_table()
{
    this->correlationTable.set_size(this->M, this->M);
    this->squaredCorrelationTable.set_size(this->M, this->M);

    for (int i = 0; i < this->M; i++) {
        for (int j = 0; j < this->M; j++) {
            this->correlationTable(i, j) = this->symbols(this->bits2symbols[i]) * conj(this->symbols(this->bits2symbols[j]));
            this->squaredCorrelationTable(i, j) = norm(this->correlationTable(i, j));
        }
    }
}

complex QAM::get_correlation(const complex& symbol1, const complex& symbol2)
{
    int index1 = this->demodulate(to_cvec(real(symbol1), imag(symbol1)))[0];
    int index2 = this->demodulate(to_cvec(real(symbol2), imag(symbol2)))[0];

    it_assert_debug(0 <= index1 && index1 < this->M, "QAM::get_correlation(): Out of range of data matrix.");
    it_assert_debug(0 <= index2 && index2 < this->M, "QAM::get_correlation(): Out of range of data matrix.");

    return this->correlationTable(index1, index2);
}

double QAM::get_squared_correlation(const complex& symbol1, const complex& symbol2)
{
    int index1 = this->demodulate(to_cvec(real(symbol1), imag(symbol1)))[0];
    int index2 = this->demodulate(to_cvec(real(symbol2), imag(symbol2)))[0];

    it_assert_debug(0 <= index1 && index1 < this->M, "QAM::get_squared_correlation(): Out of range of data matrix.");
    it_assert_debug(0 <= index2 && index2 < this->M, "QAM::get_squared_correlation(): Out of range of data matrix.");

    return this->squaredCorrelationTable(index1, index2);
}

ivec QAM::demodulate_real(const cvec& symbols)
{
    ivec indexes(symbols.size());
    int closest;

    it_assert_debug(this->setup_done, "QAM::demodulate_real(): Modulator not ready.");

    for (int i = 0; i < symbols.size(); i++) {
        closest = round_i((real(symbols[i] * this->scaling_factor) + (this->L - 1)) / 2.0);
        if (closest < 0) {
            closest = 0;
        } else if (closest > (L - 1)) {
            closest = (L - 1);
        }

        indexes[i] = closest;
    }
    return indexes;
}

cvec QAM::get_adjacent_symbols(const int index)
{
    complex symbol = this->modulate(to_ivec(index))[0];
    cvec adjacentSymbols;

    for (int i = 0; i < this->symbols.size(); i++) {
        int distance = floor_i(abs(this->symbols[i] - symbol) * this->scaling_factor);

        if (distance == 2 && distance != 0) {
            adjacentSymbols.ins(adjacentSymbols.size(), this->symbols[i]);
        }
    }
    return adjacentSymbols;
}

#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include "../../Library/Histogram/Histogram.hpp"

using std::ofstream;
using std::abs;
using std::real;
using std::pow;
using std::norm;
using std::endl;
using std::round;
using itpp::vec;
using itpp::ivec;
using itpp::cvec;
using itpp::dB;
using itpp::energy;


Histogram::Histogram(int size): DATA_SIZE(size)
{
    it_assert_debug(0 < size, "Histogram::Histogram(): Invalid number of data size.");

    this->clear();
}

Histogram::Histogram(const vec& range): DATA_SIZE(range.size())
{
    it_assert_debug(0 < range.size(), "Histogram::Histogram(): Invalid number of data size.");

    this->set_range(range);
    this->clear();
}

Histogram::Histogram(double min, double max, int size): DATA_SIZE(size)
{
    it_assert_debug(0 < size, "Histogram::Histogram(): Invalid number of data size.");

    this->set_range(min, max);
    this->clear();
}

Histogram::~Histogram()
{
}

void Histogram::clear()
{
    this->counts.set_size(DATA_SIZE);
    this->counts.zeros();
}

void Histogram::set_range(double min, double max)
{
    it_assert_debug(min < max, "Histogram::Histogram(): Invalid range designated.");

    this->min = min;
    this->max = max;
    this->stride = (this->max - this->min) / (DATA_SIZE - 1);
    this->isValid = true;
}

void Histogram::set_range(const itpp::vec& range)
{
    this->set_range(range[0], range[range.size() - 1]);
}

int Histogram::get_index(double value)
{
    return round((value - this->min) / this->stride);
}

void Histogram::increment(int index)
{
    if (index < 0) {
        this->underCount++;
    } else if (index < DATA_SIZE) {
        this->counts[index]++;
    }
    this->totalCount++;
}

void Histogram::read_number(const ivec& numbers)
{
    it_assert_debug(this->isValid == true, "Histogram::read_number(): Grapgh range is not set.");

    for (int i = 0, vectorSize = numbers.size(); i < vectorSize; i++) {
        this->increment(this->get_index(numbers[i]));
    }
}

void Histogram::read_value(const itpp::vec& values)
{
    it_assert_debug(this->isValid == true, "Histogram::read_number(): Grapgh range is not set.");

    for (int i = 0, vectorSize = values.size(); i < vectorSize; i++) {
        this->increment(this->get_index(values[i]));
    }
}

void Histogram::read_amplitude(const cvec& signals)
{
    it_assert_debug(this->isValid == true, "Histogram::read_amplitude(): Grapgh range is not set.");

    for (int i = 0, vectorSize = signals.size(); i < vectorSize; i++) {
        double value = abs(signals[i]);

        this->increment(this->get_index(value));
    }
}

void Histogram::read_phase(const cvec& signals)
{
    it_assert_debug(this->isValid == true, "Histogram::read_phase(): Grapgh range is not set.");

    for (int i = 0, vectorSize = signals.size(); i < vectorSize; i++) {
        double value = arg(signals[i]);

        this->increment(this->get_index(value));
    }
}


void Histogram::read_power(const cvec& signals)
{
    it_assert_debug(this->isValid == true, "Histogram::read_power(): Grapgh range is not set.");

    for (int i = 0, vectorSize = signals.size(); i < vectorSize; i++) {
        double value = norm(signals[i]);

        this->increment(this->get_index(value));
    }
}

void Histogram::read_normalized_power(const cvec& signals)
{
    it_assert_debug(this->isValid == true, "Histogram::read_normalized_power(): Grapgh range is not set.");

    double average = energy(signals) / signals.size();

    for (int i = 0, vectorSize = signals.size(); i < vectorSize; i++) {
        double value = dB(norm(signals[i]) / average);

        this->increment(this->get_index(value));
    }
}

vec Histogram::get_discrete_histogram()
{
    it_assert_debug(this->isValid == true, "Histogram::get_discrete_histogram(): Grapgh range is not set.");

    vec hist(DATA_SIZE);
    for (int i = 0; i < DATA_SIZE; i++) {
        hist[i] = this->counts[i] / (double)this->totalCount;
    }
    return hist;
}

vec Histogram::get_continuous_histogram()
{
    it_assert_debug(this->isValid == true, "Histogram::print_continuous_histogram(): Grapgh range is not set.");

    vec hist(DATA_SIZE);
    for (int i = 0; i < DATA_SIZE; i++) {
        hist[i] = this->counts[i] / (double)this->totalCount / this->stride;
    }
    return hist;
}

void Histogram::print_discrete_histogram(ofstream& file)
{
    it_assert_debug(this->isValid == true, "Histogram::print_discrete_histogram(): Grapgh range is not set.");

    vec hist = this->get_discrete_histogram();
    for (int i = 0; i < DATA_SIZE; i++) {
        file << this->min + this->stride * i << " " << hist[i] << endl;
    }
}

void Histogram::print_continuous_histogram(ofstream& file)
{
    it_assert_debug(this->isValid == true, "Histogram::print_continuous_histogram(): Grapgh range is not set.");

    vec hist = this->get_continuous_histogram();
    for (int i = 0; i < DATA_SIZE; i++) {
        file << this->min + this->stride * i << " " << hist[i] << endl;
    }
}

void Histogram::print_CCDF(ofstream& file)
{
    it_assert_debug(this->isValid == true, "Histogram::print_CCDF(): Grapgh range is not set.");

    double sum = this->underCount;
    for (int i = 0; i < DATA_SIZE; i++) {
        sum += this->counts[i];
        file << this->min + this->stride * i << " " << (this->totalCount - sum) / this->totalCount << endl;
    }
}

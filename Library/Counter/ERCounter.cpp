#include "../../Library/Counter/ERCounter.hpp"

ERCounter::ERCounter()
{
}

ERCounter::~ERCounter()
{
}

void ERCounter::clear()
{
    this->errors = 0;
    this->corrects = 0;
}

int ERCounter::get_errors() const
{
    return this->errors;
}

int ERCounter::get_corrects() const
{
    return this->corrects;
}

double ERCounter::get_errorrate() const
{
    if (this->corrects + this->errors > 0) {
        return this->errors / (double)(this->corrects + this->errors);
    }
    return 0.0;
}

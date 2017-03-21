#ifndef AWGNChannel_hpp
#define AWGNChannel_hpp

#include <itpp/itbase.h>
#include <itpp/base/random.h>

class AWGNChannel
{
    private:
        double sigma = 0.0;
        itpp::Complex_Normal_RNG generator;

    public:
        AWGNChannel();
        ~AWGNChannel();
        void set_noise(double);
        double get_noise() const;
        itpp::cvec operator()(const itpp::cvec&);

};

#endif

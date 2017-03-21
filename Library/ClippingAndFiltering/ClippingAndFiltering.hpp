#ifndef ClippingAndFiltering_hpp
#define ClippingAndFiltering_hpp

#include <itpp/itbase.h>

enum class CLIP_MODEL {SoftLimiter};

class ClippingAndFiltering
{
    private:
        CLIP_MODEL clipModel;
        double clippingRatio;
        int oversamplingFactor;

        itpp::cvec clip(const itpp::cvec&);
        itpp::cvec clip_by_soft_limitter(const itpp::cvec&);

    public:
        ClippingAndFiltering(double, int oversamplingFactor = 8);
        ~ClippingAndFiltering();
        double calculate_attenuation_factor();
        itpp::cvec clip_and_filter(const itpp::cvec&);
        itpp::cvec normalize_attenuation(const itpp::cvec&);

};

#endif

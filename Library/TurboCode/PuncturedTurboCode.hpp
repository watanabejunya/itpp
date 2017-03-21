#ifndef PuncturedTurboCode_hpp
#define PuncturedTurboCode_hpp

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include "../../Library/TurboCode/TurboCode.hpp"


class PuncturedTurboCode: public TurboCode
{
    protected:
        int puncturePeriod;
        long puncturedCodeLength;
        itpp::bmat punctureMatrix;
        bool tailBitUse1, tailBitUse2;

    public:
        PuncturedTurboCode();
        virtual ~PuncturedTurboCode();
        void set_puncture_matrix(const itpp::bmat&, bool = true, bool = true);
        itpp::bmat get_puncture_matrix();
        int get_puncture_period();
        int get_punctured_length();
        int get_tail_length();
        virtual itpp::bvec encode(const itpp::bvec&);
        virtual void encode(const itpp::bvec&, itpp::bvec&);
        virtual itpp::bvec decode(const itpp::bvec&);
        virtual void decode(const itpp::bvec&, itpp::bvec&);
        virtual itpp::bvec decode(const itpp::vec&);
        virtual void decode(const itpp::vec&, itpp::bvec&);
};

#endif

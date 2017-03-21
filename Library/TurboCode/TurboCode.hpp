#ifndef TurboCode_hpp
#define TurboCode_hpp

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/comm/rec_syst_conv_code.h>
#include "../../Library/Interleaver/Interleaver.hpp"


class TurboCode: itpp::Channel_Code
{
    protected:
        itpp::Rec_Syst_Conv_Code encoder1, encoder2;
        int constraintLength;
        Interleaver interleaver;
        int codeLength, informationLength;
        double rate;
        int iterations;
        std::string decodingMethod;

    public:
        TurboCode();
        virtual ~TurboCode();
        void set_parameters(itpp::ivec, itpp::ivec, int, const itpp::ivec&, int, const std::string& = "LOGMAX");
        void set_interleaver(const itpp::ivec&);
        void set_decoding_method(const std::string& method);
        int get_code_length() const;
        int get_infomation_length() const;
        virtual double get_rate() const;
        virtual itpp::bvec encode(const itpp::bvec&);
        virtual void encode(const itpp::bvec&, itpp::bvec&);
        virtual itpp::bvec decode(const itpp::bvec&);
        virtual void decode(const itpp::bvec&, itpp::bvec&);
        virtual itpp::bvec decode(const itpp::vec&);
        virtual void decode(const itpp::vec&, itpp::bvec&);
        void block_encode(const itpp::bvec&, itpp::bmat&, itpp::bmat&);
        void block_decode(const itpp::mat&, const itpp::mat&, itpp::bvec&);
};

#endif

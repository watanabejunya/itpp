#ifndef RecursiveConvolutionalCode_hpp
#define RecursiveConvolutionalCode_hpp

#include <itpp/itbase.h>
#include <itpp/comm/convcode.h>
#include <itpp/comm/llr.h>

class RecursiveConvolutionalCode
{
    private:
        std::string methodType;
        bool isTerminated;
        itpp::bvec backwardPolynomial, forwardPolynomial;
        int constraintLength;
        int memorySize;
        int outputBlockLength;
        int numberOfStates;
        double rate;
        itpp::imat stateTable, outputTable;
        itpp::LLR_calc_unit LLRCalculator;

        int calculate_next_state(const int, const int);
        int calculate_output(const int, const int);

    public:
        RecursiveConvolutionalCode();
        ~RecursiveConvolutionalCode();
        void set_generator_polynomials(const itpp::ivec&, const int);
        void set_method_type(const std::string&);
        void set_termination_type(const bool);
        virtual double get_rate() const;
        virtual itpp::bvec encode(const itpp::bvec&);
        virtual void encode(const itpp::bvec&, itpp::bvec&);
        void encode_tail(const itpp::bvec&, itpp::bvec&);
        void encode_trunc(const itpp::bvec&, itpp::bvec&);
        // virtual void decode(const itpp::bvec&, itpp::bvec&);
        // virtual itpp::bvec decode(const itpp::vec&);
        // virtual void decode(const itpp::vec&, itpp::bvec&);
        // void decode(const itpp::vec&, const itpp::vec&, const itpp::vec&, itpp::vec&);
        // void decode_table(const itpp::vec&, itpp::bvec&);
        // void decode_map(const itpp::vec&, itpp::bvec&);
        // void decode_log_map(const itpp::vec&, itpp::bvec&);
        // void decode_log_max(const itpp::vec&, itpp::bvec&);
        //
        // void map_decode(const itpp::vec&, const mat&, const itpp::vec&, itpp::vec&, bool = false);
        // void log_decode(const itpp::vec&, const itpp::vec&, const itpp::vec&, itpp::vec&);
        // void table_decode(const itpp::QLLRvec&, const itpp::QLLRvec&, const itpp::QLLRvec&, itpp::QLLRvec&);

};

#endif

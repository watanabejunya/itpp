#ifndef TrellisCode_hpp
#define TrellisCode_hpp

#include <itpp/itbase.h>

enum class TrellisCodeType {Generator, InverseSyndrome, ParityChecker, CompoundDecoder};

class TrellisCode
{
    private:
        TrellisCodeType codeType;
        bool isValid = false;
        bool isReady = false;
        int memorySize;
        int numberOfStates;
        int numberOfTransits;
        int inputBlockLength;
        int outputBlockLength;
        itpp::imat stateTable;
        itpp::imat outputTable;

    public:
        itpp::imat pathMemories;
        itpp::cmat outputSymbols;
        itpp::Mat<itpp::cvec> autocorrelation;
        itpp::mat nodeMetrics;
        itpp::mat forwardMetrics;
        itpp::mat backwardMetrics;

        TrellisCode();
        ~TrellisCode();
        void set_type(const TrellisCodeType);
        void init_trellis(const int);
        int get_memory_size();
        int get_number_of_states();
        int get_number_of_transits();
        int get_input_block_length();
        int get_output_block_length();
        itpp::ivec get_previous_states(const int);
        itpp::ivec get_next_states(const int);
        int get_next_state(const int, const int);
        int get_input(const int, const int);
        int get_output(const int, const int);
        itpp::bvec encode(const itpp::bvec&, int = 0);
        int trace_back_state(const int, const int);
        itpp::cvec get_symbols(int, int);
        int find_optimum_state(const int);
};

#endif

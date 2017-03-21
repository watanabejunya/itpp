#ifndef SelectedMapping_hpp
#define SelectedMapping_hpp

#include <itpp/itbase.h>

class SelectedMapping
{
    private:
        const int M;                        // # of candidate sequenes.
        const int N;                        // # of subcarriers.
        const int J;                        // Oversampling factor.
        itpp::ivec interleaver;             // Interleaver.

    public:
        SelectedMapping(int, int, int);                 // Constructor.
        ~SelectedMapping();                             // Deconstructor.
        itpp::cvec select(itpp::cvec);      // Select sequene with the lowest PAPR.
        itpp::cvec detect(itpp::cvec);      // Detect original information.
};

#endif

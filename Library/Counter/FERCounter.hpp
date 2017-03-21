#ifndef FERCounter_hpp
#define FERCounter_hpp

#include <itpp/itbase.h>
#include "../../Library/Counter/ERCounter.hpp"

class FERCounter: public ERCounter
{
    private:
        int blockSize;                                                          // Number of block sizes.

    public:
        FERCounter(int);                                                        // Constructor.
        virtual ~FERCounter();                                                  // Deconstructor.
        void count(const itpp::bvec&, const itpp::bvec&);                       // Count symbol errors.
};

#endif

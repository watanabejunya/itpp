#ifndef BERCounter_hpp
#define BERCounter_hpp

#include <itpp/itbase.h>
#include "../../Library/Counter/ERCounter.hpp"

class BERCounter: public ERCounter
{
    public:
        BERCounter();                                                           // Constructor.
        virtual ~BERCounter();                                                  // Deconstructor.
        void count(const itpp::bvec&, const itpp::bvec&);                       // Count symbol errors.
};

#endif

#ifndef RandomInterleaver_hpp
#define RandomInterleaver_hpp

#include "../../Library/Interleaver/Interleaver.hpp"

class RandomInterleaver: public Interleaver
{
    protected:
        void make_new_interleaver();

    public:
        RandomInterleaver(int);
        virtual ~RandomInterleaver();
        static itpp::ivec make_interleaver(const int);
};


#endif

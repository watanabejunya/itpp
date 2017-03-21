#ifndef SRandomInterleaver_hpp
#define SRandomInterleaver_hpp

#include "../../Library/Interleaver/Interleaver.hpp"

class SRandomInterleaver: public Interleaver
{
    private:
        int s;

    protected:
        void make_new_interleaver();

    public:
        SRandomInterleaver(const int, const int = 0);
        virtual ~SRandomInterleaver();
        static itpp::ivec make_interleaver(const int, const int = 0);
};


#endif

#ifndef BlockInterleaver_hpp
#define BlockInterleaver_hpp

#include "../../Library/Interleaver/Interleaver.hpp"

class BlockInterleaver: public Interleaver
{
    private:
        int rows;
        int columns;

    protected:
        void make_new_interleaver();

    public:
        BlockInterleaver(const int, const int);
        virtual ~BlockInterleaver();
        static itpp::ivec make_interleaver(const int, const int);
};


#endif

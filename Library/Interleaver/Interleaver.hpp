#ifndef Interleaver_hpp
#define Interleaver_hpp

#include <itpp/itbase.h>

class Interleaver
{
    protected:
        int size;
        itpp::ivec interleaver;

        virtual void make_new_interleaver();

    public:
        Interleaver();
        virtual ~Interleaver();
        int get_size() const;
        itpp::ivec get_interleaver() const;
        void set_interleaver(const itpp::ivec&);
        itpp::bvec interleave(const itpp::bvec&);
        itpp::vec interleave(const itpp::vec&);
        itpp::bvec deinterleave(const itpp::bvec&);
        itpp::vec deinterleave(const itpp::vec&);
};


#endif

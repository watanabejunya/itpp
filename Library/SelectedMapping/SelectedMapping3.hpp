#ifndef SelectedMapping3_hpp
#define SelectedMapping3_hpp

#include <itpp/itbase.h>
#include <itpp/comm/ofdm.h>

class SelectedMapping3
{
    private:
        const int M;                                    // # of candidate sequenes.
        const int N;                                    // # of subcarriers.
        const int J;                                    // Oversampling factor.
        itpp::cvec g1, g2, g3, g4;                      // Generator vector.
        itpp::cvec a1, a2, a3, a4;                      // Base vectors.
        itpp::cvec u1, u2, u3, u4;                      // Complex factors of unit magnitude.
        int m1, m2, m3, m4;                             // # of right syclic shifts.
        std::complex<double> p1, p2, p3, p4;            // Seeder elements.
        itpp::OFDM OFDM;                                      // OFDM modulator.

        void setGenerators();                           // Set up for generators initially.
        void setUnitFactors();                          // Set up for unit factores initially.
        void setCyclicShifts();                         // Set random cyclic shifts for seeding candidates.
        void setSeeders();                              // Set random elements for seeding candidates.
        void makeOriginVector(const itpp::cvec&);       // Make original base vector.
        void makeBaseVectors(const itpp::cvec&);        // Make base vectors.
        itpp::cvec makeCandidates();                    // Make candidate sequenes with configured parameters.

    public:
        SelectedMapping3(int m, int n, int j);                      // Constructor.
        ~SelectedMapping3();                                        // Deconstructor.
        itpp::cvec select(const itpp::cvec&);           // Select sequene with the lowest PAPR.
        itpp::cvec detect(const itpp::cvec&);           // Detect original information.
};

#endif

#include <iostream>
#include <complex>
#include <limits>
#include <itpp/itbase.h>
#include <itpp/signal/transforms.h>
#include "../../Library/Common/helper.hpp"
#include "../../Library/SelectedMapping/SelectedMapping3.hpp"

using itpp::cvec;
using itpp::randi;
using array = itpp::Array<itpp::cvec>;
using helper::mod;
using helper::rotate_left;
using helper::circular_convolute;
using helper::oversample;
using helper::ifft_normalized;
using helper::assume_PAPR;

const std::complex<double> I = std::complex<double>(0, 1);
const double INFTY = std::numeric_limits<double>::infinity();

SelectedMapping3::SelectedMapping3(int m, int n, int j): M(m), N(n), J(j)
{
    // Init settings.
    this->OFDM.set_parameters(n, 0, j);
    this->setGenerators();
    this->setUnitFactors();
}

SelectedMapping3::~SelectedMapping3()
{
}

void SelectedMapping3::setGenerators()
{
    this->g1.set_size(J * N);
    this->g2.set_size(J * N);
    this->g3.set_size(J * N);
    this->g4.set_size(J * N);
    this->g1.zeros();
    this->g2.zeros();
    this->g3.zeros();
    this->g4.zeros();

    this->g1[0] = 1;
    this->g1[J * N / 4] = 1;
    this->g1[2 * J * N / 4] = 1;
    this->g1[3 * J * N / 4] = 1;

    this->g2[0] = 1;
    this->g2[J * N / 4] = I;
    this->g2[2 * J * N / 4] = 1;
    this->g2[3 * J * N / 4] = -I;

    this->g3[0] = 1;
    this->g3[J * N / 4] = -1;
    this->g3[2 * J * N / 4] = 1;
    this->g3[3 * J * N / 4] = -1;

    this->g4[0] = 1;
    this->g4[J * N / 4] = -I;
    this->g4[2 * J * N / 4] = -1;
    this->g4[3 * J * N / 4] = I;
}

void SelectedMapping3::setUnitFactors()
{
    this->u1.set_size(4);
    this->u2.set_size(4);
    this->u3.set_size(4);
    this->u4.set_size(4);

    this->u1[0] = 1;
    this->u1[1] = 1;
    this->u1[2] = 1;
    this->u1[3] = 1;

    this->u2[0] = 1;
    this->u2[1] = I;
    this->u2[2] = -1;
    this->u2[3] = -I;

    this->u3[0] = 1;
    this->u3[1] = -1;
    this->u3[2] = 1;
    this->u3[3] = -1;

    this->u4[0] = 1;
    this->u4[1] = -I;
    this->u4[2] = -1;
    this->u4[3] = I;
}

void SelectedMapping3::setCyclicShifts()
{
    this->m1 = 0;
    this->m2 = randi(0, J * N / 4 - 1);
    this->m3 = randi(0, J * N / 4 - 1);
    this->m4 = randi(0, J * N / 4 - 1);
}

void SelectedMapping3::setSeeders()
{
    cvec seeds(4);
    seeds[0] = 1;
    seeds[1] = -1;
    seeds[2] = I;
    seeds[3] = -I;

    this->p1 = 1;
    this->p2 = seeds[randi(0, 3)];
    this->p3 = seeds[randi(0, 3)];
    this->p4 = seeds[randi(0, 3)];
}

void SelectedMapping3::makeBaseVectors(const cvec& signals)
{
    this->a1.set_size(J * N);
    this->a2.set_size(J * N);
    this->a3.set_size(J * N);
    this->a4.set_size(J * N);

    this->a1.set_subvector(0, circular_convolute(signals, this->g1, J * N / 4));
    this->a1.set_subvector(J * N / 4, this->u1[1] * this->a1.left(J * N / 4));
    this->a1.set_subvector(J * N * 2 / 4, this->u1[2] * this->a1.left(J * N / 4));
    this->a1.set_subvector(J * N * 3 / 4, this->u1[3] * this->a1.left(J * N / 4));

    this->a2.set_subvector(0, circular_convolute(signals, this->g2, J * N / 4));
    this->a2.set_subvector(J * N / 4, this->u2[1] * this->a2.left(J * N / 4));
    this->a2.set_subvector(J * N * 2 / 4, this->u2[2] * this->a2.left(J * N / 4));
    this->a2.set_subvector(J * N * 3 / 4, this->u2[3] * this->a2.left(J * N / 4));

    this->a3.set_subvector(0, circular_convolute(signals, this->g3, J * N / 4));
    this->a3.set_subvector(J * N / 4, this->u3[1] * this->a3.left(J * N / 4));
    this->a3.set_subvector(J * N * 2 / 4, this->u3[2] * this->a3.left(J * N / 4));
    this->a3.set_subvector(J * N * 3 / 4, this->u3[3] * this->a3.left(J * N / 4));

    this->a4.set_subvector(0, circular_convolute(signals, this->g4, J * N / 4));
    this->a4.set_subvector(J * N / 4, this->u4[1] * this->a4.left(J * N / 4));
    this->a4.set_subvector(J * N * 2 / 4, this->u4[2] * this->a4.left(J * N / 4));
    this->a4.set_subvector(J * N * 3 / 4, this->u4[3] * this->a4.left(J * N / 4));
}

cvec SelectedMapping3::makeCandidates()
{
    return this->a1 + this->p2 * rotate_left(this->a2, this->m2)
        + this->p3 * rotate_left(this->a3, this->m3) + this->p4 * rotate_left(this->a4, this->m4) ;
}

cvec SelectedMapping3::select(const cvec& symbols)
{
    cvec signals = this->OFDM.modulate(symbols);

    this->makeBaseVectors(signals);

    array candidates(M);
    for (int i = 0; i < M; i++) {
        // Prepare for making candidates.
        this->setCyclicShifts();
        this->setSeeders();

        candidates(i) = this->makeCandidates();
    }

    // Select symbols with the lowest PAPR.
    double lowestPAPR = INFTY;
    int needle;
    for (int i = 0; i < M; i++) {
        double PAPR = assume_PAPR(candidates(i));

        // Update the lowest one.
        if(PAPR < lowestPAPR) {
            lowestPAPR = PAPR;
            needle = i;
        }
    }

    return this->OFDM.demodulate(candidates(needle));
}

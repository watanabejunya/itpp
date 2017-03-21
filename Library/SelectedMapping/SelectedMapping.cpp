#include <iostream>
#include <complex>
#include <itpp/itbase.h>
#include "../../Library/SelectedMapping/SelectedMapping.hpp"

using itpp::cvec;

SelectedMapping::SelectedMapping(int m, int n, int j): M(m), N(n), J(j)
{
}

SelectedMapping::~SelectedMapping()
{
}

cvec SelectedMapping::select(cvec)
{

}

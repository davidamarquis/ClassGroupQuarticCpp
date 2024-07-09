#ifndef THESIS_PREDEFINEDIDEALS_HPP
#define THESIS_PREDEFINEDIDEALS_HPP

#include "src/Quartic/QuarticField.hpp"
#include "src/Quartic/QuarticIdeal.hpp"
#include "src/IndexCalculus/Quartic/FactorBase.hpp"

inline const char* small_nf_str = "-2, 0, 0, 0, 1";
inline const char* non_monogenic_nf_str = "1, 1, 2, -1, 1";

void SmallQuarticNf(QuarticField &nf);
void SmallQuarticNfFactorBase(FactorBase &fb, QuarticField &nf);

void IdlInSmallQuarticNf(QuarticIdeal &idl, QuarticField &nf);
void IdlWithBasisZeroInSmallQuarticNf(QuarticIdeal &idl, QuarticField &nf);

void NonMonogenicNf(QuarticField &nf);

void IdlInNonMonogenicNf(QuarticIdeal &idl, QuarticField &nf);
void IdlsInNonMonogenicNf(QuarticIdeal &idl0, QuarticIdeal &idl1, QuarticField &nf);

void IdlsWithBasisZeroInNonMonogenicNf(QuarticIdeal &idl0, QuarticIdeal &idl1, QuarticField &nf);
void IdlWithBasisZeroInNonMonogenicNf(QuarticIdeal &idl, QuarticField &nf);
#endif  // THESIS_PREDEFINEDIDEALS_HPP
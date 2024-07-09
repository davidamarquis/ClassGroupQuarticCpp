#include "PredefinedIdeals.hpp"
#include "src/IndexCalculus/Quartic/FactorBase.hpp"

void SmallQuarticNf(QuarticField &nf) {
  // nf=Q[X]/(X^4-2)
  QuarticFieldFromCParams(nf, small_nf_str, "1", "0, 1", "0, 0, 1",
                          "0, 0, 0, 1", "2048", "1, 0, 0, 0\n0, 1, 0, 0\n0, 0, 1, 0\n0, 0, 0, 1", "1");
}

void IdlInSmallQuarticNf(QuarticIdeal &idl, QuarticField &nf) {
  SmallQuarticNf(nf);
  // idl=(7, Y-2)
  QuarticIdealFromCParams(idl, nf, 7, "-2, 1","7, 0, 0, 0\n5, 1, 0, 0\n3, 0, 1, 0\n6, 0, 0, 1", 1, 1);
}

void IdlWithBasisZeroInSmallQuarticNf(QuarticIdeal &idl, QuarticField &nf) {
  SmallQuarticNf(nf);
  // idl=(7, Y-2)
  QuarticIdealFromCParams(idl, nf, 7, "-2, 1", nullptr, 1, 1);
}

void SmallQuarticNfFactorBase(FactorBase &fb, QuarticField &nf) {
  SmallQuarticNf(nf);
  uint32_t smoothness_bound = 5000;
  fb = FactorBaseFromJson("DataForTests/SmallQuarticNfFacBase.json", nf, smoothness_bound);
}

void NonMonogenicNf(QuarticField &nf) {
  // nf=Q[X]/(X^4 - X^3 + 2*X^2 + X + 1)
  QuarticFieldFromCParams(nf,non_monogenic_nf_str, "1", "-1/2, 1, -1, 1/2", "0, 1", "1/2, 0, 0, 1/2", "225", "1, 0, 0, 0\n0, 0, 1, 0\n-1, -1, 1, 1\n-1, 0, 0, 2", "1");
}

void IdlInNonMonogenicNf(QuarticIdeal &idl, QuarticField &nf) {
  NonMonogenicNf(nf);
  // idl=(19, Y-2)
  QuarticIdealFromCParams(idl, nf, 19, "-2, 1", "19, 0, 0, 0\n8, 1, 0, 0\n17, 0, 1, 0\n5, 0, 0, 1", 1, 1);
}

void IdlsInNonMonogenicNf(QuarticIdeal &idl0, QuarticIdeal &idl1, QuarticField &nf) {
  IdlInNonMonogenicNf(idl0, nf);
  QuarticIdealFromCParams(idl1, nf, 31, "2, 1", "31, 0, 0, 0\n26, 1, 0, 0\n2, 0, 1, 0\n19, 0, 0, 1", 1, 1);
}

void IdlWithBasisZeroInNonMonogenicNf(QuarticIdeal &idl, QuarticField &nf) {
  NonMonogenicNf(nf);
  // idl=(7, Y-2)
  QuarticIdealFromCParams(idl, nf, 19, "-2, 1", nullptr, 1, 1);
}

void IdlsWithBasisZeroInNonMonogenicNf(QuarticIdeal &idl0, QuarticIdeal &idl1, QuarticField &nf) {
  IdlWithBasisZeroInNonMonogenicNf(idl0, nf);
  QuarticIdealFromCParams(idl1, nf, 31, "2, 1", nullptr, 1, 1);
}


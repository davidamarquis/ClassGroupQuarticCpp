#include <src/Quartic/PredefinedIdeals.hpp>
#include <src/Quartic/QuarticField.hpp>
#include <src/IndexCalculus/Quartic/FactorBase.hpp>

#include "tests/catch.hpp"
using namespace std;

void test_TrivialRelations() {
  const char* expected_mat_str="4, 0, 0, 0, 0, 0, 0\n0, 1, 1, 0, 0, 0, 0\n0, 0, 0, 1, 0, 0, 0\n0, 0, 0, 0, 1, 1, 1";
  fmpz_mat_t expected_mat;
  fmpz_mat_init(expected_mat,4,7); // degree=4, number of primes=7
  mat_str_to_fmpz_mat(expected_mat, expected_mat_str);

  const char* fb_path = "DataForTests/SmallQuarticNfFacBase.json";
  QuarticField nf = QuarticField(small_nf_str);
  SmallQuarticNf(nf);
  FactorBase fb = FactorBaseFromJson(fb_path, nf, 7);
  fmpz_mat_t mat;
  uint32_t nrows = fb.rat_primes.size();
  uint32_t ncols = fb.prime_idls.size();
  fmpz_mat_init(mat, nrows, ncols);

  fb.TrivialRelations(mat);
  REQUIRE(fmpz_mat_equal(mat, expected_mat));
}


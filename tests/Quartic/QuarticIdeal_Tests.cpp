#include <flint/flint.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include <cstdio>
#include <iostream>
#include <src/Flint.hpp>
#include <src/IndexCalculus/Quartic/ObjStringConverters.hpp>
#include <src/Quartic/PredefinedIdeals.hpp>
#include <src/Quartic/QuarticField.hpp>
#include <src/Quartic/QuarticIdeal.hpp>

#include "tests/catch.hpp"

using namespace std;

void ValidateIdlIntegralBasisCoordinates(QuarticIdeal &idl, QuarticField &nf, const char* expected_nf_bas_coords_mat_num_str, int expected_nf_bas_coords_mat_dnm, const char* expected_idl_bas_coords_mat_num_str, int expected_idl_bas_coords_mat_dnm) {
  fmpz_mat_t bas_coords_z;
  fmpq_mat_t bas_coords;
  fmpq_mat_t prod;
  nf_elem_t elem;
  fmpz_mat_init(bas_coords_z,4,4);
  fmpq_mat_init(bas_coords,4,4);
  fmpq_mat_init(prod,4,4);

  nf_elem_init(elem, nf.antic_nf);
  fmpz_t dnms[4];
  for (int i=0;i<4;i++) {
    fmpz_init(dnms[i]);
  }

//  field tests
  for (int i=0;i<4;i++) {
    nf.GetBasisElem(i, elem);
    nf_elem_get_fmpz_mat_row(bas_coords_z, i, dnms[i], elem, nf.antic_nf);
  }
  fmpq_mat_set_fmpz_mat(bas_coords, bas_coords_z);
  fmpq_mat_divide_rows(bas_coords, dnms);

  fmpq_mat_t expected_nf_bas_coords_mat;
  fmpq_mat_init(expected_nf_bas_coords_mat,4,4);
  num_mat_str_and_den_int_to_fmpq_mat(expected_nf_bas_coords_mat, expected_nf_bas_coords_mat_dnm, expected_nf_bas_coords_mat_num_str);
  CHECK(fmpq_mat_equal(expected_nf_bas_coords_mat, bas_coords));
  fmpq_mat_clear(expected_nf_bas_coords_mat);

  fmpq_mat_t id4;
  fmpq_mat_init(id4,4,4);
  fmpq_mat_one(id4);
  fmpq_mat_mul(prod, bas_coords, nf.int_basis_mat_inv);
  CHECK(fmpq_mat_equal(prod, id4));
  fmpq_mat_clear(id4);

  // ideal tests
  for (int i=0;i<4;i++) {
    idl.GetBasisElem(i, elem);
    nf_elem_get_fmpz_mat_row(bas_coords_z, i, dnms[i], elem, nf.antic_nf);
  }
  fmpq_mat_set_fmpz_mat(bas_coords, bas_coords_z);
  fmpq_mat_divide_rows(bas_coords, dnms);

  fmpq_mat_t expected_idl_bas_coords_mat;
  fmpq_mat_init(expected_idl_bas_coords_mat,4,4);
  num_mat_str_and_den_int_to_fmpq_mat(expected_idl_bas_coords_mat, expected_idl_bas_coords_mat_dnm, expected_idl_bas_coords_mat_num_str);
  CHECK(fmpq_mat_equal(bas_coords, expected_idl_bas_coords_mat));
  fmpq_mat_clear(expected_idl_bas_coords_mat);

  fmpq_mat_mul(prod,  bas_coords, nf.int_basis_mat_inv);
  CHECK(fmpq_mat_equal_fmpz_mat(prod, idl.basis_mat));

  fmpz_mat_clear(bas_coords_z);
  fmpq_mat_clear(bas_coords);
  fmpq_mat_clear(prod);
  nf_elem_clear(elem, nf.antic_nf);
  for (int i=0;i<4;i++) {
    fmpz_clear(dnms[i]);
  }
}

TEST_CASE("Ideal - in two calls to RandomElement the elements returned are not the same") {
  QuarticField nf = QuarticField(small_nf_str);
  flint_rand_t state;
  flint_randinit(state);
  SECTION("Basis Nonzero") {
    QuarticIdeal idl = QuarticIdeal(nf);
    IdlInSmallQuarticNf(idl, nf);

    nf_elem_t rand_elem0;
    nf_elem_init(rand_elem0, nf.antic_nf);
    nf_elem_t rand_elem1;
    nf_elem_init(rand_elem1, nf.antic_nf);

    idl.RandomElem(rand_elem0, state);
    idl.RandomElem(rand_elem1, state);
    CHECK(!nf_elem_equal(rand_elem0, rand_elem1, nf.antic_nf));
    nf_elem_clear(rand_elem0, nf.antic_nf);
    nf_elem_clear(rand_elem1, nf.antic_nf);
  }

  SECTION("Basis Matrix is Zero") {
    QuarticIdeal idlb = QuarticIdeal(nf);
    IdlWithBasisZeroInSmallQuarticNf(idlb, nf);
    CHECK(fmpz_mat_is_zero(idlb.basis_mat));
    CHECK(nf_elem_is_zero(idlb.b0, nf.antic_nf));

    nf_elem_t rand_elem0b;
    nf_elem_init(rand_elem0b, nf.antic_nf);
    nf_elem_t rand_elem1b;
    nf_elem_init(rand_elem1b, nf.antic_nf);

    idlb.RandomElemWithGens(rand_elem0b, state);
    idlb.RandomElemWithGens(rand_elem1b, state);
    CHECK(!nf_elem_equal(rand_elem0b, rand_elem1b, nf.antic_nf));
    nf_elem_clear(rand_elem0b, nf.antic_nf);
    nf_elem_clear(rand_elem1b, nf.antic_nf);
  }

  flint_randclear(state);
}

TEST_CASE("Ideal - tau0", "[Quartic]") {
  QuarticField nf = QuarticField("1, 1, 2, -1, 1");
  QuarticIdeal idl = QuarticIdeal(nf);
  IdlInSmallQuarticNf(idl, nf);
  const char* uniformizer_str = "-2, 1";
  nf_elem_t uniformizer; // a non-random uniformizer, just for this test
  nf_elem_init(uniformizer, nf.antic_nf);
  poly_str_to_nf_elem(uniformizer, uniformizer_str, nf.antic_nf);

  fmpq_mat_t uniformizer_mul_mat;
  fmpq_mat_init(uniformizer_mul_mat,4,4);

  CHECK(nf_elem_is_one(idl.uniformizer, nf.antic_nf));
  CHECK(nf_elem_is_one(idl.tau0, nf.antic_nf));
  flint_rand_t state;
  flint_randinit(state);
  nf_elem_set(idl.uniformizer, uniformizer, nf.antic_nf);
  idl.InitTau0(uniformizer_mul_mat, state);
  assert(nf_elem_equal(idl.uniformizer, uniformizer, nf.antic_nf));

  const char* expected_mat_str = "-2, 1, 0, 0\n0, -2, 1, 0\n0, 0, -2, 1\n2, 0, 0, -2";
  fmpq_mat_t expected_uniformizer_mul_mat;
  fmpq_mat_init(expected_uniformizer_mul_mat,4,4);
  num_mat_str_and_den_int_to_fmpq_mat(expected_uniformizer_mul_mat, 1, expected_mat_str);

  CHECK(fmpq_mat_equal(uniformizer_mul_mat, expected_uniformizer_mul_mat));
  nf_elem_t red_tau0;
  nf_elem_init(red_tau0, nf.antic_nf);
  nf_elem_set_si(red_tau0, 23, nf.antic_nf); // set the initial value of tau0 to something selected randomly
  nf_elem_mul(red_tau0, idl.tau0, idl.uniformizer, nf.antic_nf);
// TODO get coeffs of numerator of tau0 (mod p)
//  nf_elem_mod_fmpz(red_tau0, red_tau0, idl.gen0, nf.antic_nf);
//
//  CHECK(nf_elem_is_zero(red_tau0, nf.antic_nf));
//  nf_elem_clear(red_tau0, nf.antic_nf);
}

TEST_CASE("Ideal - tau0 non-monogenic", "[Quartic]") {
  QuarticField nf = QuarticField(non_monogenic_nf_str);
  QuarticIdeal idl = QuarticIdeal(nf);
  IdlInNonMonogenicNf(idl, nf);
  const char* uniformizer_str = "360, 6, 0, 1";
  nf_elem_t uniformizer; // a non-random uniformizer, just for this test
  nf_elem_init(uniformizer, nf.antic_nf);
  poly_str_to_nf_elem(uniformizer, uniformizer_str, nf.antic_nf);

  fmpq_mat_t uniformizer_mul_mat;
  fmpq_mat_init(uniformizer_mul_mat,4,4);

  CHECK(nf_elem_is_one(idl.uniformizer, nf.antic_nf));
  CHECK(nf_elem_is_one(idl.tau0, nf.antic_nf));
  flint_rand_t state;
  flint_randinit(state);
  nf_elem_set(idl.uniformizer, uniformizer, nf.antic_nf);
  idl.InitTau0(uniformizer_mul_mat, state);
  assert(nf_elem_equal(idl.uniformizer, uniformizer, nf.antic_nf));

  const char* expected_mat_str = "359, 0, 6, 2\n0, 359, -4, -6\n-6, -4, 363, 6\n2, 6, -6, 357";
  fmpq_mat_t expected_uniformizer_mul_mat;
  fmpq_mat_init(expected_uniformizer_mul_mat,4,4);
  num_mat_str_and_den_int_to_fmpq_mat(expected_uniformizer_mul_mat, 1, expected_mat_str);

  CHECK(fmpq_mat_equal(uniformizer_mul_mat, expected_uniformizer_mul_mat));
  nf_elem_t red_tau0;
  nf_elem_init(red_tau0, nf.antic_nf);
  nf_elem_set_si(red_tau0, 23, nf.antic_nf); // set the initial value of tau0 to something selected randomly
  nf_elem_mul(red_tau0, idl.tau0, idl.uniformizer, nf.antic_nf);
// TODO get coeffs of numerator of tau0 (mod p)
//
//  CHECK(nf_elem_is_zero(red_tau0, nf.antic_nf));
//  nf_elem_clear(red_tau0, nf.antic_nf);
}

void ValuationTestsForNfAndIdeals(QuarticField &nf, QuarticIdeal &idl0, QuarticIdeal &idl1) {

  fmpq_mat_t bas_mat;
  fmpq_mat_init(bas_mat,4,4);
  flint_rand_t rand;
  flint_randinit(rand);
  idl0.InitTau0(bas_mat, rand);
  idl1.InitTau0(bas_mat, rand);
  fmpq_mat_clear(bas_mat);
  flint_randclear(rand);

  SECTION("Valuations equal to 0 and to 1") {
    nf_elem_t elem0;
    nf_elem_init(elem0, nf.antic_nf);
    nf_elem_set(elem0, idl0.gen1, nf.antic_nf);
    int val0i0 = idl0.Valuation(elem0);
    int val0i1 = idl1.Valuation(elem0);
    CHECK(val0i0 == 1);
    CHECK(val0i1 == 0);
    nf_elem_clear(elem0, nf.antic_nf);
  }

  SECTION("Valuation is larger than one") {
    nf_elem_t elem1;
    nf_elem_init(elem1, nf.antic_nf);
    nf_elem_pow(elem1, idl0.gen1, 2, nf.antic_nf);
    int val1 = idl0.Valuation(elem1);
    CHECK(val1 == 2);
    nf_elem_clear(elem1, nf.antic_nf);
  }

  SECTION("Principal ideal generated by an element factors into two prime ideals") {
    nf_elem_t elem2;
    nf_elem_init(elem2, nf.antic_nf);
    nf_elem_pow(elem2, idl0.gen1, 2, nf.antic_nf);
    nf_elem_mul(elem2, elem2, idl1.gen1, nf.antic_nf);
    int val2i0 = idl0.Valuation(elem2);
    int val2i1 = idl1.Valuation(elem2);
    CHECK(val2i0 == 2);
    CHECK(val2i1 == 1);
    nf_elem_clear(elem2, nf.antic_nf);
  }
}

TEST_CASE("Ideal - Valuation Non-monogenic Field", "[Quartic]") {
  QuarticField nf = QuarticField(non_monogenic_nf_str);
  QuarticIdeal idl0 = QuarticIdeal(nf);
  QuarticIdeal idl1 = QuarticIdeal(nf);
  IdlsInNonMonogenicNf(idl0, idl1, nf);
  ValuationTestsForNfAndIdeals(nf, idl0, idl1);
}

TEST_CASE("Ideal - Valuation Non-monogenic Field Basis Mat Equal To 0", "[Quartic]") {
  QuarticField nf = QuarticField(non_monogenic_nf_str);
  QuarticIdeal idl0 = QuarticIdeal(nf);
  QuarticIdeal idl1 = QuarticIdeal(nf);
  IdlsWithBasisZeroInNonMonogenicNf(idl0, idl1, nf);
  ValuationTestsForNfAndIdeals(nf, idl0, idl1);
}

TEST_CASE("Ideal - Integral Basis Coords", "[Quartic]") {
  SECTION("Small field") {
    QuarticField nf = QuarticField("1, 1, 2, -1, 1");
    QuarticIdeal idl = QuarticIdeal(nf);
    IdlInSmallQuarticNf(idl, nf);
    ValidateIdlIntegralBasisCoordinates(idl, nf, "1, 0, 0, 0\n0, 1, 0, 0\n0, 0, 1, 0\n0, 0, 0, 1",1, "7, 0, 0, 0\n5, 1, 0, 0\n3, 0, 1, 0\n6, 0, 0, 1", 1); // good
  }
  SECTION("Non-monogenic field") {
    QuarticField nf2 = QuarticField("1, 1, 2, -1, 1");
    QuarticIdeal idl2 = QuarticIdeal(nf2);
    IdlInNonMonogenicNf(idl2, nf2);
    ValidateIdlIntegralBasisCoordinates(
        idl2, nf2, "2, 0, 0, 0\n-1, 2, -2, 1\n0, 2, 0, 0\n1, 0, 0, 1", 2,
        "38, 0, 0, 0\n15, 2, -2, 1\n34, 2, 0, 0\n11, 0, 0, 1", 2);
  }
}

TEST_CASE("Ideal - basis", "[Quartic]") {
  QuarticField nf = QuarticField("1, 1, 2, -1, 1");
  QuarticIdeal idl = QuarticIdeal(nf);
  IdlInSmallQuarticNf(idl, nf);

  fmpq_poly_t expected_b0_pol;
  fmpq_poly_init(expected_b0_pol);
  fmpq_poly_set_str(expected_b0_pol, "1  7");
  fmpq_poly_t expected_b1_pol;
  fmpq_poly_init(expected_b1_pol);
  fmpq_poly_set_str(expected_b1_pol, "2  5 1");
  fmpq_poly_t expected_b2_pol;
  fmpq_poly_init(expected_b2_pol);
  fmpq_poly_set_str(expected_b2_pol, "3  3 0 1");
  fmpq_poly_t expected_b3_pol;
  fmpq_poly_init(expected_b3_pol);
  fmpq_poly_set_str(expected_b3_pol, "4  6 0 0 1");

  CHECK(fmpq_poly_equal(NF_ELEM(idl.b0), expected_b0_pol));
  CHECK(fmpq_poly_equal(NF_ELEM(idl.b1), expected_b1_pol));
  CHECK(fmpq_poly_equal(NF_ELEM(idl.b2), expected_b2_pol));
  CHECK(fmpq_poly_equal(NF_ELEM(idl.b3), expected_b3_pol));
  fmpq_poly_clear(expected_b0_pol);
  fmpq_poly_clear(expected_b1_pol);
  fmpq_poly_clear(expected_b2_pol);
  fmpq_poly_clear(expected_b3_pol);
}

#include <cstdio>
#include <iostream>
#include <src/Quartic/PredefinedIdeals.hpp>
#include <src/Quartic/QuarticField.hpp>
#include <src/Quartic/QuarticIdeal.hpp>

#include "tests/catch.hpp"
using namespace std;

TEST_CASE("Quartic Field - ElemToOkCoordinates", "[Quartic]") {
  QuarticField nf = QuarticField(non_monogenic_nf_str);
  QuarticIdeal idl0 = QuarticIdeal(nf);
  QuarticIdeal idl1 = QuarticIdeal(nf);
  IdlsInNonMonogenicNf(idl0, idl1, nf);
  fmpz_mat_t coord_mat0;
  fmpz_mat_init(coord_mat0, 1, 4);
  bool ret0 = nf.ElemToOkCoordinates(coord_mat0, idl0.gen1);
  CHECK(ret0);

  fmpz_mat_t coord_mat1;
  fmpz_mat_init(coord_mat1, 1, 4);
  bool ret1 = nf.ElemToOkCoordinates(coord_mat1, idl1.gen1);
  CHECK(ret1);
  fmpz_mat_clear(coord_mat1);

  fmpz_mat_t coord_mat2;
  fmpz_mat_init(coord_mat2, 1, 4);
  nf_elem_t elem2;
  nf_elem_init(elem2, nf.antic_nf);
  nf_elem_set_si(elem2, 1, nf.antic_nf);
  nf_elem_div(elem2, elem2, idl0.gen1, nf.antic_nf);
  bool ret2 = nf.ElemToOkCoordinates(coord_mat0, elem2);
  CHECK(!ret2);
  nf_elem_clear(elem2, nf.antic_nf);
  fmpz_mat_clear(coord_mat2);

  fmpz_mat_t coord_mat3;
  fmpz_mat_init(coord_mat3, 1, 4);
  nf_elem_t elem3;
  nf_elem_init(elem3, nf.antic_nf);
  nf_elem_set_si(elem3, 0, nf.antic_nf);
  bool ret3 = nf.ElemToOkCoordinates(coord_mat3, elem3);
  CHECK(ret3);
  CHECK(fmpz_mat_is_zero(coord_mat3));
  fmpz_mat_clear(coord_mat3);
  nf_elem_clear(elem3, nf.antic_nf);

  bool ret4;
  fmpz_mat_t coord_mat4;
  fmpz_mat_init(coord_mat4, 1, 4);
  nf_elem_t elem4;
  nf_elem_init(elem4, nf.antic_nf);
  for (int cc=0; cc<4; cc++) {
    nf.GetBasisElem(cc, elem4);
    ret4 = nf.ElemToOkCoordinates(coord_mat4, elem4);
    CHECK(ret4);
    CHECK(fmpz_is_one(fmpz_mat_entry(coord_mat4, 0, cc)));
  }
  fmpz_mat_clear(coord_mat4);
  nf_elem_clear(elem4, nf.antic_nf);
}


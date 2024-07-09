#include "QuarticField.hpp"
#include <cassert>

#include "include/rapidjson/document.h"
#include "include/rapidjson/istreamwrapper.h"
#include "include/rapidjson/ostreamwrapper.h"
#include "include/rapidjson/stringbuffer.h"
#include "include/rapidjson/writer.h"
#include "src/IndexCalculus/Quartic/ObjStringConverters.hpp"
#include "src/Quartic/QuarticIdeal.hpp"
#include "src/Flint.hpp"
#include "src/IndexCalculus/Quartic/IoUtils.hpp"

using namespace rapidjson;
using namespace std;

bool QuarticField::operator==(const QuarticField &other) {
  return fmpz_poly_equal(def_poly, other.def_poly) and nf_elem_equal(w0, other.w0, antic_nf) and nf_elem_equal(w1, other.w1, antic_nf) and nf_elem_equal(w2, other.w2, antic_nf) and nf_elem_equal(w3, other.w3, antic_nf) and fmpz_equal(Delta, other.Delta);
}

void QuarticField::InitValueMemberVariables() {
  // must be called after antic_nf is initialized
  nf_elem_init(w0, antic_nf);
  nf_elem_init(w1, antic_nf);
  nf_elem_init(w2, antic_nf);
  nf_elem_init(w3, antic_nf);
  fmpz_init(Delta);
  fmpz_poly_init(def_poly);
  fmpq_mat_init(int_basis_mat_inv,4,4);
}

void QuarticField::assign(const QuarticField &nf) {
  // looks like reinitialization is required for the nf_elems. It's simpler to just reinitialize all member variables
  Clear();

  fmpq_poly_t def_poly_q;
  fmpq_poly_init(def_poly_q);
  fmpq_poly_set_fmpz_poly(def_poly_q, nf.def_poly);
  nf_init(antic_nf, def_poly_q);
  fmpq_poly_clear(def_poly_q);

  InitValueMemberVariables();
  nf_elem_set(w0, nf.w0, antic_nf);
  nf_elem_set(w1, nf.w1, antic_nf);
  nf_elem_set(w2, nf.w2, antic_nf);
  nf_elem_set(w3, nf.w3, antic_nf);
  fmpz_set(Delta, nf.Delta);

  fmpz_poly_set(def_poly, nf.def_poly);
  fmpq_mat_set(int_basis_mat_inv, nf.int_basis_mat_inv);
}

QuarticField::QuarticField(const QuarticField&& nf) {
  if (!fmpz_poly_equal(def_poly, nf.def_poly)) {
    throw std::invalid_argument("Assignment is only supported for fields with the same defining polynomial");
  }
  assign(nf);
}

QuarticField::QuarticField(const QuarticField& nf) {
  if (!fmpz_poly_equal(def_poly, nf.def_poly)) {
    throw std::invalid_argument("Assignment is only supported for fields with the same defining polynomial");
  }
  assign(nf);
}

QuarticField::QuarticField(const char* poly_str) {
  /*
   * Initializes the field with a valid antic_nf and all other member variables default initialized
   */
  fmpq_poly_t def_poly;
  fmpq_poly_init(def_poly);
  poly_str_to_fmpq_poly(def_poly, poly_str);

  if (fmpq_poly_degree(def_poly) != 4) {
    throw std::invalid_argument("Defining poly must have degree 4");
  }
  nf_init(antic_nf, def_poly); // valid antic_nf
  fmpq_poly_clear(def_poly);

  InitValueMemberVariables(); // default initialization
}

QuarticField::QuarticField(fmpq_poly_t def_poly_in, fmpq_poly_t w0_in, fmpq_poly_t w1_in, fmpq_poly_t w2_in, fmpq_poly_t w3_in, fmpz_t Delta_in, fmpq_mat_t int_basis_mat_inv_in, vector<int64_t> poly_coeffs_in) {
  if (fmpq_poly_degree(def_poly_in) != 4) {
    throw std::invalid_argument("Defining poly must have degree 4");
  }
  fmpq_mat_init(int_basis_mat_inv,4,4);
  fmpq_mat_set(int_basis_mat_inv, int_basis_mat_inv_in);
  copy(poly_coeffs_in.begin(), poly_coeffs_in.end(), back_inserter(poly_coeffs));

  nf_init(antic_nf, def_poly_in);
  fmpz_init(Delta);
  fmpz_set(Delta, Delta_in);

  fmpz_t den;
  fmpz_init(den);
  fmpq_poly_get_denominator(den, def_poly_in);
  if (!fmpz_is_one(den)) {
    throw std::invalid_argument("Defining poly must have coeffs in ZZ");
  }
  if (!fmpq_poly_is_monic(def_poly_in)) {
    throw std::invalid_argument("Defining poly must be monic");
  }
  fmpz_poly_init(def_poly);
  fmpq_poly_get_numerator(def_poly, def_poly_in);
#ifndef NDEBUG
  fmpz_t disc;
  fmpz_init(disc);
  fmpz_poly_discriminant(disc, def_poly);
  if (!fmpz_divisible(disc, Delta)) {
    throw std::invalid_argument("Delta must divide the defining poly's discriminant");
  }
#endif

  nf_elem_init(w0, antic_nf);
  nf_elem_set_fmpq_poly(w0, w0_in, antic_nf);
  nf_elem_init(w1, antic_nf);
  nf_elem_set_fmpq_poly(w1, w1_in, antic_nf);
  nf_elem_init(w2, antic_nf);
  nf_elem_set_fmpq_poly(w2, w2_in, antic_nf);
  nf_elem_init(w3, antic_nf);
  nf_elem_set_fmpq_poly(w3, w3_in, antic_nf);
}

QuarticField& QuarticField::operator=(const QuarticField &other) {
  if (this != &other) {
    assign(other);
  }
  return *this;
}

QuarticField& QuarticField::operator=(const QuarticField &&other) {
  if (this != &other) {
    assign(other);
  }
  return *this;
}

void QuarticField::GetBasisElem(int i, nf_elem_t bas_elem) {
  if (i==0) {
    nf_elem_set(bas_elem, w0, antic_nf);
  } else if (i==1) {
    nf_elem_set(bas_elem, w1, antic_nf);
  } else if (i==2) {
    nf_elem_set(bas_elem, w2, antic_nf);
  } else if (i==3) {
    nf_elem_set(bas_elem, w3, antic_nf);
  } else {
    assert(false);
  }
}

void QuarticField::Clear() {
  nf_elem_clear(w0, antic_nf);
  nf_elem_clear(w1, antic_nf);
  nf_elem_clear(w2, antic_nf);
  nf_elem_clear(w3, antic_nf);
  fmpz_poly_clear(def_poly);
  fmpq_mat_clear(int_basis_mat_inv);
  fmpz_clear(Delta);

  nf_clear(antic_nf);
}

QuarticField::~QuarticField() {
  printf("Destroy field ");
  this->print();

  Clear();
}

void QuarticField::fprint(FILE *fp) {
  fprintf(fp, "Q[X]/(");
  fmpz_poly_fprint_pretty(fp, def_poly,"X");
  fprintf(fp, ")\n");
}

void QuarticField::print() {
  fprint(stdout);
}

void QuarticField::OkCoordinateMatRowToElem(nf_elem_t new_elem, const fmpz_mat_t coord_mat, const slong rr) {
  /*
   * Takes a mat whose rows are the coordinates of elems of Ok on the integral basis and returns the element corresponding to row rr
   */
  nf_elem_t tmp_elem;
  nf_elem_t bb;
  nf_elem_init(tmp_elem, antic_nf);
  nf_elem_init(bb, antic_nf);
  nf_elem_set_si(new_elem, 0, antic_nf);
  for (slong cc=0; cc<fmpz_mat_ncols(coord_mat); cc++) {
    nf_elem_set_fmpz(tmp_elem, fmpz_mat_entry(coord_mat, rr, cc), antic_nf);
    GetBasisElem(cc, bb);
    nf_elem_mul(tmp_elem, tmp_elem, bb, antic_nf);
    nf_elem_add(new_elem, new_elem, tmp_elem, antic_nf);
  }
  nf_elem_clear(tmp_elem, antic_nf);
  nf_elem_clear(bb, antic_nf);
}

bool QuarticField::ElemToOkCoordinates(fmpz_mat_t coord_mat, const nf_elem_t elem) {
/*
 * returns false if elem not in Ok and leaves coord_mat unchanged
 * returns true if the elem is in Ok. This can be used to test if an element is on the integral basis
 */
  if (fmpz_mat_nrows(coord_mat)!=1 or fmpz_mat_ncols(coord_mat)!=4) {
    throw std::invalid_argument("coord_mat must have dimensions 1x4");
  }
  fmpz_mat_t coords;
  fmpz_mat_init(coords,1,4);
  fmpq_mat_t coords_q;
  fmpq_mat_init(coords_q,1,4);
  fmpz_t elem_dnm;
  fmpz_init(elem_dnm);

  nf_elem_get_fmpz_mat_row(coords, 0, elem_dnm, elem, antic_nf);
  fmpq_mat_set_fmpz_mat(coords_q, coords);  // "cast" coords to QQ entries
  fmpq_mat_scalar_div_fmpz(coords_q, coords_q, elem_dnm);

  fmpq_mat_mul(coords_q, coords_q, int_basis_mat_inv);

  fmpz_t num;
  fmpz_init(num);
  fmpz_t dnm;
  fmpz_init(dnm);

  bool ret = true;
  for (int cc=0;cc<4;cc++) {
    fmpq_canonicalise(fmpq_mat_entry(coords_q, 0, cc));
  }
  ret = fmpq_mat_is_integral(coords_q);
  if (ret) {
    fmpq_mat_get_fmpz_mat(coord_mat, coords_q);
  }

  fmpz_clear(dnm);
  fmpz_clear(elem_dnm);
  fmpz_mat_clear(coords);
  fmpq_mat_clear(coords_q);
  return ret;
}

void QuarticField::RatPrimeToQuarticIdeals(vector<QuarticIdeal> &idls, const ulong prime) {
  /*
   * Results are undefined if prime is not a prime
   */
  nmod_poly_t poly;
  nmod_poly_init(poly, prime);
  fmpz_poly_get_nmod_poly(poly, def_poly);
  if (fmpz_divisible_si(Delta, prime)) {
    throw std::invalid_argument("prime cannot divide Delta");
  }
  if (!nmod_poly_is_squarefree(poly)) {
    throw std::invalid_argument("Defining polynomial must be squarefree (mod p).");
  }
  if (nmod_poly_is_irreducible(poly)) {
    throw std::invalid_argument("Defining polynomial is irreducible (mod p). P is inert in the number field.");
  }

  nmod_poly_factor_t facs;
  nmod_poly_factor_init(facs);
  nmod_poly_factor(facs, poly);

  fmpz_t gen0;
  fmpz_init_set_ui(gen0, prime);
  fmpq_poly_t elem_pol;
  fmpq_poly_init(elem_pol);
  fmpz_mat_t zero;
  fmpz_mat_init(zero,4,4);
  fmpz_mat_zero(zero);
  QuarticIdeal idl = QuarticIdeal(*this);

  for(int i=0;i<facs->num;i++) {
    nmod_poly_print(&(facs->p)[i]);
    fmpq_poly_set_nmod_poly(elem_pol, &(facs->p)[i]);
    idl = QuarticIdeal(gen0, elem_pol, zero, true, nmod_poly_degree(&(facs->p)[i]), (facs->exp)[i], *this);
    idls.push_back(idl);
  }
  nmod_poly_factor_clear(facs);
  fmpz_clear(gen0);
  fmpq_poly_clear(elem_pol);
  fmpz_mat_clear(zero);
  nmod_poly_clear(poly);
}

void QuarticFieldFromCParams(QuarticField &nf, const char* def_poly_str, const char* w0_str, const char* w1_str, const char* w2_str, const char* w3_str, const char* Delta_str, const char* int_basis_mat_inv_num_str, const char* int_basis_mat_inv_dnm_str) {
  fmpq_poly_t ff;
  fmpq_poly_init(ff);
  poly_str_to_fmpq_poly(ff, def_poly_str);

  // set integral basis for this monogenic field
  fmpq_poly_t w0;
  fmpq_poly_init(w0);
  poly_str_to_fmpq_poly(w0, w0_str);

  fmpq_poly_t w1;
  fmpq_poly_init(w1);
  poly_str_to_fmpq_poly(w1, w1_str);

  fmpq_poly_t w2;
  fmpq_poly_init(w2);
  poly_str_to_fmpq_poly(w2, w2_str);

  fmpq_poly_t w3;
  fmpq_poly_init(w3);
  poly_str_to_fmpq_poly(w3, w3_str);

  fmpz_t Delta;
  fmpz_init(Delta);
  fmpz_set_str(Delta, Delta_str, 10);

  fmpz_mat_t int_basis_mat_inv_num;
  fmpz_t int_basis_mat_inv_dnm;
  fmpq_mat_t int_basis_mat_inv;
  fmpz_mat_init(int_basis_mat_inv_num, 4, 4);
  fmpz_init(int_basis_mat_inv_dnm);
  fmpq_mat_init(int_basis_mat_inv, 4, 4);
  mat_str_to_fmpz_mat(int_basis_mat_inv_num, int_basis_mat_inv_num_str);
  fmpz_set_str(int_basis_mat_inv_dnm, int_basis_mat_inv_dnm_str, 10);
  fmpq_mat_set_fmpz_mat_div_fmpz(int_basis_mat_inv, int_basis_mat_inv_num, int_basis_mat_inv_dnm);

  vector<string> coeffs_strs {};
  Split(",", string(def_poly_str), coeffs_strs);
  vector<int64_t> coeffs {};
  for (auto s: coeffs_strs) {
    coeffs.push_back(stol(s));
  }

  nf = QuarticField(ff,w0,w1,w2,w3,Delta, int_basis_mat_inv, coeffs);

  fmpq_poly_clear(ff);
  fmpq_poly_clear(w0);
  fmpq_poly_clear(w1);
  fmpq_poly_clear(w2);
  fmpq_poly_clear(w3);
  fmpz_clear(Delta);
  fmpz_mat_clear(int_basis_mat_inv_num);
  fmpz_clear(int_basis_mat_inv_dnm);
  fmpq_mat_clear(int_basis_mat_inv);
}

QuarticField QuarticFieldFromJson(const string nf_json_path) {
  Document nf_doc{};
  LoadDoc(nf_doc, nf_json_path);

  GenericArray pari_int_basis = nf_doc["pari_int_bas"]["basis"].GetArray();
  QuarticField nf = QuarticField(nf_doc["def_poly"].GetString());
  const char *poly_str;
  poly_str = nf_doc["def_poly"].GetString();

  QuarticFieldFromCParams(
      nf, poly_str, pari_int_basis[0].GetString(),
      pari_int_basis[1].GetString(), pari_int_basis[2].GetString(),
      pari_int_basis[3].GetString(), nf_doc["Delta"].GetString(),
      nf_doc["int_basis_mat_inv_num"].GetString(),
      nf_doc["int_basis_mat_inv_dnm"].GetString());
  return nf;
}


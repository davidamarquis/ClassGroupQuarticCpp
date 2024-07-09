#ifndef THESIS_QUARTICFIELD_HPP
#define THESIS_QUARTICFIELD_HPP

#include <iostream>
#include <stdexcept>
#include <flint/flint.h>
#include <flint/arith.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq_poly.h>
#include <antic/nf.h>
#include <antic/nf_elem.h>
#include <src/IndexCalculus/Quartic/ObjStringConverters.hpp>
#include <src/Quartic/QuarticIdeal.hpp>

class QuarticIdeal;

class QuarticField {
  public:
  int dg = 4;
  nf_elem_t w0;
  nf_elem_t w1;
  nf_elem_t w2;
  nf_elem_t w3;
  fmpz_t Delta; //absolute value of the field's discrim
  nf_t antic_nf;
  fmpz_poly_t def_poly; // defining polynomials are required to be integral
  fmpq_mat_t int_basis_mat_inv;
  vector<int64_t> poly_coeffs;

  bool operator==(const QuarticField &other);

  void InitValueMemberVariables();

  void assign(const QuarticField &nf);

  QuarticField() = delete;
  QuarticField(const QuarticField&& nf);

  QuarticField(const QuarticField& nf);

  QuarticField(const char* poly_str); // minimal initialization
  QuarticField(fmpq_poly_t def_poly_in, fmpq_poly_t w0_in, fmpq_poly_t w1_in, fmpq_poly_t w2_in, fmpq_poly_t w3_in, fmpz_t Delta_in, fmpq_mat_t int_basis_mat_inv_in, vector<int64_t> poly_coeffs_in);

  QuarticField& operator=(const QuarticField &other);

  QuarticField& operator=(const QuarticField &&other);

  void GetBasisElem(int i, nf_elem_t bas_elem);
  void OkCoordinateMatRowToElem(nf_elem_t new_elem, const fmpz_mat_t coord_mat, const slong rr);
  bool ElemToOkCoordinates(fmpz_mat_t coord_mat, const nf_elem_t elem);

  void Clear();
  ~QuarticField();

  void fprint(FILE *fp);

  void print();

  void RatPrimeToQuarticIdeals(vector<QuarticIdeal> &idls, const ulong prime);
};

void QuarticFieldFromCParams(QuarticField &nf, const char* def_poly_str, const char* w0_str, const char* w1_str, const char* w2_str, const char* w3_str, const char* Delta_str, const char* int_basis_mat_inv_num_str, const char* int_basis_mat_inv_dnm_str);
QuarticField QuarticFieldFromJson(const string nf_json_path);

#endif  // THESIS_QUARTICFIELD_HPP
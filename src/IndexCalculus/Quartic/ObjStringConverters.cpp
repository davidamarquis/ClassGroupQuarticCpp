#include <algorithm>

#include <flint/arith.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include "ObjStringConverters.hpp"
#include <src/IndexCalculus/Quartic/IoUtils.hpp>

using namespace std;

void poly_str_to_coeffs(vector<string> &coeffs, string &poly) {
  /*
   * Initializes an fmpq_poly corresponding to the input string.
   * The format of poly_c_str is a sequence of integers representing the coefficients 0 through d and separated by commas.
   */
}

void poly_str_to_nmod_poly(nmod_poly_t &poly, const char *poly_c_str) {
  /*
   * Initializes an fmpq_poly corresponding to the input string.
   * The format of poly_c_str is a sequence of integers representing the coefficients 0 through d and separated by commas.
   */
  string poly_str {poly_c_str};
  auto dg = std::count(poly_str.begin(), poly_str.end(), ',');
  auto poly_size = dg + 1;

  poly_str.erase(remove(poly_str.begin(), poly_str.end(), ','), poly_str.end());

  // append length of the poly
  poly_str = to_string(poly_size) + "  " + poly_str;
  // reduce poly coeffs
  nmod_poly_set_str(poly, poly_str.c_str());
}

void poly_str_to_fmpq_poly(fmpq_poly_t poly, const char *poly_c_str) {
  /*
   * Initializes an fmpq_poly corresponding to the input string.
   * The format of poly_c_str is a sequence of integers representing the coefficients 0 through d and separated by commas.
   */
  string poly_str {poly_c_str};
  auto dg = std::count(poly_str.begin(), poly_str.end(), ',');
  auto poly_size = dg + 1;

  poly_str.erase(remove(poly_str.begin(), poly_str.end(), ','), poly_str.end());

  // append length of the poly
  poly_str = to_string(poly_size) + "  " + poly_str;
  fmpq_poly_set_str(poly, poly_str.c_str());
}

void poly_str_to_nf_elem(nf_elem_t elem, const char *poly_c_str, const nf_t nf) {
  /*
   * Sets elem to the nf_elem defined by poly_c_str.
   * The format of poly_c_str is the same as poly_str_to_fmpq_poly().
   */
  fmpq_poly_t poly;
  fmpq_poly_init(poly);
  poly_str_to_fmpq_poly(poly, poly_c_str);
  nf_elem_set_fmpq_poly(elem, poly, nf);
  fmpq_poly_clear(poly);
}

void flint_mat_str_to_fmpz_mat(fmpz_mat_t mat, const string str) {
  FILE* tmpf = std::tmpfile();
  fputs(str.c_str(), tmpf);
  rewind(tmpf);
  fmpz_mat_fread(tmpf, mat);
}

void mat_str_to_flint_mat_str(string &flint_mat_str, const char *mat_str) {
  if (mat_str == nullptr) {
    throw runtime_error("Matrix string is null");
  }

  string delimiter = "\n";
  string mat_strpp = string(mat_str);
  vector<string> mat_row_strs {};
  Split(delimiter, mat_strpp, mat_row_strs);

  int num_rows = mat_row_strs.size();
  string::difference_type num = count(mat_row_strs[0].begin(), mat_row_strs[0].end(), ',');
  int num_cols = num+1;
  flint_mat_str = to_string(num_rows);
  flint_mat_str += " ";
  flint_mat_str += to_string(num_cols);
  flint_mat_str += "  ";
  string entry_str = mat_row_strs[0];
  for (int rr = 1; rr < num_rows; rr++) {
    entry_str += " ";
    entry_str += mat_row_strs[rr];
  }
  flint_mat_str += entry_str;
  flint_mat_str.erase(remove(flint_mat_str.begin(), flint_mat_str.end(), ','), flint_mat_str.end());
}

void mat_str_to_fmpz_mat(fmpz_mat_t mat, const char *mat_str) {
  if (mat_str == nullptr) {
    throw runtime_error("Matrix string is null");
  }
  string flint_mat_str = "";
  mat_str_to_flint_mat_str(flint_mat_str, mat_str);
  flint_mat_str_to_fmpz_mat(mat, flint_mat_str);
}

void num_mat_str_and_den_int_to_fmpq_mat(fmpq_mat_t mat, const int den_int, const char *num_mat_str) {
  if (num_mat_str == nullptr) {
    throw runtime_error("Numerator matrix string is null");
  }
  fmpz_mat_t num;
  fmpz_mat_init(num,4,4);
  mat_str_to_fmpz_mat(num, num_mat_str);

  fmpz_t den;
  fmpz_init_set_si(den, den_int);

  fmpq_mat_set_fmpz_mat_div_fmpz(mat, num, den);
}

#include <src/Flint.hpp>
#include <stdexcept>

void nf_elem_norm_abs(fmpz_t nrm_abs, nf_elem_t elem, nf_t nf) {
  fmpq_t nrm;
  fmpq_init(nrm);
  nf_elem_norm(nrm, elem, nf);
  fmpz_abs(nrm_abs, &(nrm->num));
  fmpq_clear(nrm);
}

slong nf_elem_norm_valn(fmpz_t prime, nf_elem_t elem, nf_t nf) {
  fmpq_t nrm;
  fmpz_t nrm_z;
  nf_elem_norm(nrm, elem, nf);

  fmpq_numerator(nrm_z, nrm);
  // If this function needs to be made faster consider replacing fmpz_remove with calling _fmpz_remove with a precomputed inverse
  slong ss = fmpz_remove(nrm_z, nrm_z, prime);

  fmpq_clear(nrm);
  fmpz_clear(nrm_z);
  return ss;
}

int fmpq_mat_equal_fmpz_mat(const fmpq_mat_t q_mat, const fmpz_mat_t z_mat) {
  fmpq_mat_t mat;
  fmpq_mat_init(mat,4,4);
  fmpq_mat_set_fmpz_mat(mat,z_mat);
  int res = fmpq_mat_equal(mat, q_mat);
  fmpq_mat_clear(mat);
  return res;
}

void fmpq_mat_divide_rows(fmpq_mat_t mat, const fmpz_t *dnms) {
  /*
   * Intended usage is for mat to be a matrix with the denominator of each row cleared into the vector dnms
   */
  slong nrows = fmpq_mat_nrows(mat);
  slong ncols = fmpq_mat_ncols(mat);
  if (nrows != 4) {
    throw std::invalid_argument("mat must have four rows");
  }

  fmpq_t entry;
  fmpq_init(entry);
  for (int rr = 0; rr < nrows; rr++) {
    for (int cc = 0; cc < ncols; cc++) {
      fmpq_div_fmpz(entry, fmpq_mat_entry(mat, rr, cc), dnms[rr]);
      fmpq_set(fmpq_mat_entry(mat, rr, cc), entry);
    }
  }
  fmpq_clear(entry);
}

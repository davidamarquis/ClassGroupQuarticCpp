#ifndef THESIS_FLINT_HPP
#define THESIS_FLINT_HPP

#include <flint/flint.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq_mat.h>
#include <antic/nf.h>
#include <antic/nf_elem.h>

int fmpq_mat_equal_fmpz_mat(const fmpq_mat_t q_mat, const fmpz_mat_t z_mat);

void fmpq_mat_divide_rows(fmpq_mat_t mat, const fmpz_t *dnms);

slong nf_elem_norm_valn(fmpz_t prime, nf_elem_t elem, nf_t nf);

void nf_elem_norm_abs(fmpz_t nrm, nf_elem_t, nf_t);
#endif
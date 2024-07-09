#ifndef THESIS_CONVERTERS_HPP
#define THESIS_CONVERTERS_HPP

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <iomanip>

#include <ctime>
#include <fstream>
#include <cstdlib>
#include <string>

#include <flint/arith.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <antic/nf.h>
#include <antic/nf_elem.h>

using namespace std;

void poly_str_to_nmod_poly(nmod_poly_t &poly, const char *poly_c_str);

void poly_str_to_fmpq_poly(fmpq_poly_t poly, const char *poly_c_str);

void poly_str_to_nf_elem(nf_elem_t elem, const char *poly_c_str, const nf_t nf);

void flint_mat_str_to_fmpz_mat(fmpz_mat_t mat, const string str);

void mat_str_to_flint_mat_str(string &flint_mat_str, const char *mat_str);

void mat_str_to_fmpz_mat(fmpz_mat_t mat, const char *mat_str);

void num_mat_str_and_den_int_to_fmpq_mat(fmpq_mat_t mat, const int den_int, const char *num_mat_str);

#endif
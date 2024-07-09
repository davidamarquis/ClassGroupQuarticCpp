#pragma once

#include <set>
#include "flint/arith.h"
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

bool row_eq(slong row_ind0, fmpz_mat_t const(&mat0), slong row_ind1,
            fmpz_mat_t const(&mat1));
void inds_of_dup_rows(const fmpz_mat_t &mat, std::set<ulong> &inds_dups, std::set<ulong> &inds_of_non_dup_rows);
void swap_rows_to_top(fmpz_mat_t &mat, slong *perm, const std::set<ulong> inds_to_swap, std::set<ulong> other_inds);
void remove_dups(fmpz_mat_t &mat, fmpz_mat_t &submat);


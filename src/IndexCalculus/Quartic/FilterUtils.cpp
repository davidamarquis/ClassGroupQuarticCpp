#include <vector>
#include <iostream>
#include <set>
#include <cstdio>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <cassert>
#include "FilterUtils.hpp"

using namespace std;

bool row_eq(slong row_ind0, fmpz_mat_t const(&mat0), slong row_ind1,
            fmpz_mat_t const(&mat1)) {
  /*
   * Returns mat0[row_ind0] == mat1[row_ind1]
   */
  slong ncols = fmpz_mat_ncols(mat0);
  fmpz_t m_0_ri0_l;
  fmpz_t m_1_ri1_l;
  bool ret_val = true;
  for (slong l=0; l<ncols; l++) {
    //TODO get rid of init here
    fmpz_init_set(m_0_ri0_l, fmpz_mat_entry(mat0, row_ind0, l));
    fmpz_init_set(m_1_ri1_l, fmpz_mat_entry(mat1, row_ind1, l));
    if (*m_0_ri0_l != *m_1_ri1_l) {
      ret_val = false;
      break;
    }
  }
  fmpz_clear(m_0_ri0_l);
  fmpz_clear(m_1_ri1_l);
  return ret_val;
}

void remove_dups(fmpz_mat_t &mat, fmpz_mat_t &submat) {
  /*
   * initializes submat be the submatrix of all the unique rows of mat
   */
  slong nrows = fmpz_mat_nrows(mat);
  slong ncols = fmpz_mat_ncols(mat);
  set<ulong> inds_dups;
  set<ulong> non_dup_rows;
  inds_of_dup_rows(mat, inds_dups, non_dup_rows);
  swap_rows_to_top(mat, NULL, inds_dups, non_dup_rows);
  slong r1 = inds_dups.size();
  fmpz_mat_window_init(submat, mat, r1, 0, nrows, ncols);
}

void swap_rows_to_top(fmpz_mat_t &mat, slong *perm, const std::set<ulong> inds_to_swap, std::set<ulong> other_inds) {
  // modifies mat and other_inds in-place. The resulting other_inds is not the same length and does not contain meaningful data
  for (ulong i : inds_to_swap) {
    if (i < inds_to_swap.size()) {
      continue;
    }
    ulong dest_ind = *other_inds.begin(); // first entry in other_inds
    fmpz_mat_swap_rows(mat, NULL, i, dest_ind);
    other_inds.erase(dest_ind);
  }
}

void inds_of_dup_rows(const fmpz_mat_t &mat, set<ulong> &inds_dups, set<ulong> &inds_of_non_dup_rows) {
  //Computes a partition [inds_of_dup_rows, inds_of_non_dup_rows] of [0,...,nrows-1] so that for every i in
  // inds_of_dup_rows there is a j in inds_of_non_dup_rows so that mat[i]=mat[j]
  //:param mat:
  //:return:
//  vector<fmpz_> visited_rows = vector<int>();
  ulong nrows = fmpz_mat_nrows(mat);
  inds_dups = set<ulong>();

  inds_of_non_dup_rows = set<ulong>();
  for(ulong i = 0; i < nrows; i++) {
    inds_of_non_dup_rows.insert(i);
  }

  for (ulong j=0; j<nrows; j++) {
    for (ulong k = 0; k < j; k++) {
      if (row_eq(j, mat, k, mat)) {
        inds_dups.insert(j);
        inds_of_non_dup_rows.erase(j);
      }
    }
  }
}

void remove_dup_timing() {
  // checks remove dups takes the right amount of time
  fmpz_mat_t mat;
  fmpz_mat_init(mat, 0, 0);
  FILE *filePointer;
  filePointer = fopen("DataForTests/flint_mat3", "r");
  if (filePointer == NULL) {
    assert(false);
  }
  fmpz_mat_fread(filePointer, mat);
  fclose(filePointer);
  fmpz_mat_t submat;

  clock_t c_start = clock();

  remove_dups(mat, submat);

  clock_t c_end = clock();
  clock_t time_in_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  cout << fixed << setprecision(2) << "CPU time used: " << time_in_ms << " ms\n";

  set<ulong> inds_dups2;
  set<ulong> non_dup_rows;
  inds_of_dup_rows(submat, inds_dups2, non_dup_rows);
  assert(inds_dups2.size()==0);
  assert(time_in_ms < 200);
  assert(time_in_ms > 100);

  fmpz_mat_clear(mat);
  fmpz_mat_clear(submat);
}


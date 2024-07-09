#include <iostream>
#include <src/IndexCalculus/Quartic/FilterUtils.hpp>

#include "tests/catch.hpp"

using namespace std;

TEST_CASE("Filter - row_eq", "[IndCalc]") {
  FILE *fp = fopen("DataForTests/flint_mat0", "r");
  if (fp == NULL) {
    REQUIRE(false);
  }
  fmpz_mat_t mat;
  fmpz_mat_init(mat, 2, 2);
  fmpz_mat_fread(fp, mat);
  slong ncols = fmpz_mat_ncols(mat);
  REQUIRE(row_eq(0, mat, 0, mat));
  REQUIRE(row_eq(1, mat, 1, mat));
  REQUIRE(row_eq(1, mat, 0, mat));
}

TEST_CASE("Filter - inds_of_dup_rows 0", "[IndCalc]") {
  FILE *fp = fopen("DataForTests/flint_mat0", "r");
  if (fp == NULL) {
    REQUIRE(false);
  }
  fmpz_mat_t test_mat;
  fmpz_mat_init(test_mat, 2, 2);
  fmpz_mat_fread(fp, test_mat);
  set<ulong> inds_dups;
  set<ulong> inds_of_non_dup_rows;
  inds_of_dup_rows(test_mat, inds_dups, inds_of_non_dup_rows);
  REQUIRE(inds_dups.size()==1);
  REQUIRE(inds_of_non_dup_rows.size()==1);
  fmpz_mat_clear(test_mat);
  fclose(fp);
}

TEST_CASE("Filter - inds_of_dup_rows 1", "[IndCalc]") {
  FILE *fp = fopen("DataForTests/flint_mat1", "r");
  if (fp == NULL) {
    REQUIRE(false);
  }
  fmpz_mat_t test_mat;
  fmpz_mat_init(test_mat, 5, 2);
  fmpz_mat_fread(fp, test_mat);
  set<ulong> inds_dups;
  set<ulong> inds_of_non_dup_rows;
  inds_of_dup_rows(test_mat, inds_dups, inds_of_non_dup_rows);
  REQUIRE(inds_dups.size()==2);
  REQUIRE(inds_of_non_dup_rows.size()==3);
  REQUIRE(std::find(inds_dups.begin(), inds_dups.end(), 1) != inds_dups.end());
  REQUIRE(std::find(inds_dups.begin(), inds_dups.end(), 4) != inds_dups.end());
  fmpz_mat_clear(test_mat);
  fclose(fp);
}

TEST_CASE("Filter - swap rows 1", "[IndCalc]") {
//  [[1,2],[1,2],[3,4],[5,6],[1,2]]
  FILE *fp = fopen("DataForTests/flint_mat1", "r");
  if (fp == NULL) {
    REQUIRE(false);
  }
  fmpz_mat_t test_mat;
  fmpz_mat_init(test_mat, 5, 2);
  fmpz_mat_fread(fp, test_mat);
  fmpz_mat_t expected_mat;
  fmpz_mat_init_set(expected_mat, test_mat);

  set<ulong> inds_dups;
  set<ulong> inds_of_non_dup_rows;
  inds_of_dup_rows(test_mat, inds_dups, inds_of_non_dup_rows);
  slong *perm;
  swap_rows_to_top(test_mat, perm, inds_dups, inds_of_non_dup_rows);
  SECTION("Inds of dups is correct") {
    REQUIRE(inds_dups.size() == 2);
    REQUIRE(inds_of_non_dup_rows.size() == 3);
    REQUIRE(std::find(inds_dups.begin(), inds_dups.end(), 1) !=
            inds_dups.end());
    REQUIRE(std::find(inds_dups.begin(), inds_dups.end(), 4) !=
            inds_dups.end());
  }
  SECTION("Mat's rows are unchanged") {
    REQUIRE(fmpz_mat_equal(test_mat, expected_mat));
  }
  fmpz_mat_clear(test_mat);
  fmpz_mat_clear(expected_mat);
  fclose(fp);
}

TEST_CASE("Filter - swap rows 2", "[IndCalc]") {
//  [[3,4],[5,6],[1,2],[1,2],[1,2]]
  FILE *fp = fopen("DataForTests/flint_mat2", "r");
  if (fp == NULL) {
    REQUIRE(false);
  }
  fmpz_mat_t test_mat;
  fmpz_mat_init(test_mat, 5, 2);
  fmpz_mat_fread(fp, test_mat);
  fclose(fp);

  FILE *fp2 = fopen("DataForTests/flint_mat2_dups_top", "r");
  if (fp2 == NULL) {
    REQUIRE(false);
  }
  fmpz_mat_t expected_mat;
  fmpz_mat_init(expected_mat, 5, 2);
  fmpz_mat_fread(fp, expected_mat);
  fclose(fp2);

  set<ulong> inds_dups;
  set<ulong> inds_of_non_dup_rows;
  inds_of_dup_rows(test_mat, inds_dups, inds_of_non_dup_rows);
  slong *perm;
  swap_rows_to_top(test_mat, perm, inds_dups, inds_of_non_dup_rows);
  SECTION("Mat's rows are unchanged") {
    REQUIRE(fmpz_mat_equal(test_mat, expected_mat));
  }
  fmpz_mat_clear(test_mat);
  fmpz_mat_clear(expected_mat);
}

TEST_CASE("Filter - swap rows, 8.2mb matrix", "[IndCalc]") {
  FILE *fp = fopen("DataForTests/flint_mat3", "r");
  if (fp == NULL) {
    REQUIRE(false);
  }
  fmpz_mat_t test_mat;
  fmpz_mat_init(test_mat, 0, 0);
  fmpz_mat_fread(fp, test_mat);
  fclose(fp);

  set<ulong> inds_dups;
  set<ulong> inds_of_non_dup_rows;
  inds_of_dup_rows(test_mat, inds_dups, inds_of_non_dup_rows);
  slong *perm;
  swap_rows_to_top(test_mat, perm, inds_dups, inds_of_non_dup_rows);
  slong ncols = fmpz_mat_ncols(test_mat);
  fmpz_mat_t window;
  fmpz_mat_window_init(window, test_mat, 0, 0, 1, ncols);
  SECTION("First row is correct") {
    REQUIRE(*fmpz_mat_entry(window, 0, 0) == 0);
    REQUIRE(*fmpz_mat_entry(window, 0, 1) == 0);
    REQUIRE(*fmpz_mat_entry(window, 0, 2) == 4);
    REQUIRE(*fmpz_mat_entry(window, 0, 3) == 0);
    REQUIRE(*fmpz_mat_entry(window, 0, 4) == 1);
    REQUIRE(*fmpz_mat_entry(window, 0, 5) == 0);
    REQUIRE(*fmpz_mat_entry(window, 0, 6) == 2);
    REQUIRE(*fmpz_mat_entry(window, 0, 7) == 1);
    REQUIRE(*fmpz_mat_entry(window, 0, 8) == 1);
    REQUIRE(*fmpz_mat_entry(window, 0, 9) == 0);
  }
  fmpz_mat_clear(test_mat);
  fmpz_mat_clear(window);
}
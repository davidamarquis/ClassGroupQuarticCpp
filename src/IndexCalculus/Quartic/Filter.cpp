#include "Filter.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstdint> // for fixed-width int types
#include <cassert>
#include <cmath> // for pow() function
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <antic/nf.h>
#include <antic/nf_elem.h>
#include <flint/arith.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>

#include "include/rapidjson/document.h"
#include "include/rapidjson/istreamwrapper.h"
#include "include/rapidjson/ostreamwrapper.h"
#include "include/rapidjson/stringbuffer.h"
#include "include/rapidjson/writer.h"

#include "src/Flint.hpp"
#include "src/Quartic/PredefinedIdeals.hpp"
#include "src/Quartic/QuarticField.hpp"
#include "src/Quartic/QuarticElem.hpp"
#include "src/Quartic/QuarticIdeal.hpp"
#include "src/IndexCalculus/Quartic/FilterUtils.hpp"
#include "src/IndexCalculus/Quartic/IoUtils.hpp"
#include "src/IndexCalculus/Quartic/FactorBase.hpp"
#include "src/IndexCalculus/Quartic/ObjStringConverters.hpp"

inline const char* idl_coeff_gcd_str = "idl_coeff_gcd";
inline const char* idl_least_int_str = "idl_least_int";
inline const char* idl_poly_gen_str = "idl_poly_gen";
inline const char* idl_fb_exp_vec_str = "idl_fb_exp_vec";

void filter(FactorBase& fb, QuarticField& nf, uint32_t num_ideals, float filtering_ratio_of_relns_to_load, int8_t lpb, int32_t smoothness_bound) {
  const char* filename = "/Users/davidmarquis/git/thesis/ANTL/tests/DataForTests/tmp/all_reln";
  std::ifstream f(filename);
  if (!f.is_open()) {
    throw std::invalid_argument(string("Could not open filename: all_reln at: ") + string(filename));
  }
  cout << "dg of nf is " << nf.dg << endl;

  std::string line;
  long total_reln_available = 0;
  uint64_t max_num_lines = pow(2, 24);  // never going to have sieved more than 2^24=16777216 relations
  int i = 0;
  vector<string> words = vector<string>();
  vector<string> words_for_time = vector<string>();
  string delimiter{" "};
  double total_las_time = 0;

  vector<uint64_t> idl_least_ints{};
  nf_elem_t qf_elem;
  nf_elem_init(qf_elem, nf.antic_nf);
  // idls have the form <idl_least_ints, idl_poly_gens>
  vector<QuarticElem> idl_poly_gens{};
  vector<QuarticElem> quartic_elems{};
  vector<uint64_t> idl_coeff_gcds{};
  int32_t idl_index = -1;

//  sieve_pr_facs.clear();
//  Split(string(","), split_res[1], sieve_pr_facs_hex);

  vector<int32_t> norm_pr_facs = vector<int32_t>();
  vector<string> sieve_pr_facs_hex = vector<string>();

  // first pass over file
  for (i = 0; i < max_num_lines; i++) {
    // idl_index needs to be compared to an unsigned int. This complicated b/c in a comparison x_signed < y_unsigned C++ makes a weird choice and converts x_signed to an unsigned value.
    if (idl_index >= (int)num_ideals) {
      throw runtime_error("idl_index too large");
    }
    getline(f, line);
    if (f.eof()) {
      break;
    }

    // TODO start with checking the first character is '#'
    if (StartsWith("# Total cpu time", line)) {
      words_for_time.clear();
      Split(string(" "), line, words_for_time);
      total_las_time += stod(words_for_time[4].c_str(), nullptr);
    } else if (StartsWith("# Total", line)) {
      words.clear();
      Split(" ", line, words);
      if (words.size() < 4 and words[3] == "reports") {
        continue;
      }
      total_reln_available += strtol(words[2].c_str(), nullptr, 10);
    }

    if (StartsWith(idl_least_int_str, line)) {
      idl_index += 1;
      vector<string> res{};
      Split(idl_least_int_str + string(":"), line, res);
      if (res.size() < 2) {
        string err_str = string("Error when splitting the ideal with index ") + to_string(idl_index) + string(" with the str ") + string(idl_least_int_str);
        throw runtime_error(err_str);
      }
      idl_least_ints.push_back(strtol(res[1].c_str(), nullptr, 10));
    } else if (StartsWith(idl_poly_gen_str, line)) {
      vector<string> res_left_bkt{};
      vector<string> res_right_bkt{};
      Split(string("["), line, res_left_bkt);
      if (res_left_bkt.size() < 2) {
        string err_str = string("Error when splitting the ideal with index ") + to_string(idl_index) + string(" with the str ") + string("'['");
        throw runtime_error(err_str);
      }
      Split(string("]"), res_left_bkt[1], res_right_bkt);
      if (res_right_bkt.size() < 2) {
        string err_str = string("Error when splitting the ideal with index ") + to_string(idl_index) + string(" with the str ") + string("']'");
        throw runtime_error(err_str);
      }
      poly_str_to_nf_elem(qf_elem, res_right_bkt[0].c_str(), nf.antic_nf);
      idl_poly_gens.emplace_back(qf_elem);
//      idl_poly_gens.push_back(qf_elem);

    } else if (StartsWith(idl_coeff_gcd_str, line)) {
      vector<string> split_res{};
      Split(idl_coeff_gcd_str + string(":"), line, split_res);
      if (split_res.size() < 2) {
        string err_str = string("Error when splitting the ideal with index ") + to_string(idl_index) + string(" with the str ") + string(idl_coeff_gcd_str);
        throw runtime_error(err_str);
      }
      idl_coeff_gcds.push_back(strtol(split_res[1].c_str(), nullptr, 10));
    }
  }

  string stage{"Filter"};
  // sanity checks
  if (total_reln_available == 0) {
    cerr << stage << "Reading every relation from all_relns. all_reln file has {} lines" << endl;
  }
  if (idl_coeff_gcds.size() == 0) {
    throw runtime_error("idl_coeff_gcds is empty");
  }
  float rat_lines_used = i;
  rat_lines_used /= max_num_lines;
  cout << "Total las time" << total_las_time << endl;
  cout << "Fraction of max_num_lines used " << std::setprecision(2) << std::fixed << rat_lines_used << endl;
  cout << "total reln available " << total_reln_available << " vs BB= " << fb.fb_bound << endl;
  cout << "num lines examined in file " << i << endl;

  cout << "first in idl_least_ints" << idl_least_ints.at(0) << endl;
  //  cout << "first in idl_poly_gens" << idl_poly_gens[0] << endl;
  cout << "first in idl_coeff_gcds" << idl_coeff_gcds.at(0) << endl;

  int32_t num_prs = fb.prime_idls.size();
  int32_t num_relns_requested = floor(filtering_ratio_of_relns_to_load * fb.prime_idls.size());
  num_relns_requested = total_reln_available;  // TODO remove
  if (num_relns_requested > total_reln_available) {
    vector<bool> initial_reln_inds = vector<bool>(total_reln_available, 1);
    num_relns_requested = total_reln_available;
    cout << "Reading every relation from all_relns" << endl;
  } else {
    // get a random 0-1 tuple that's guaranteed to be at least as long as the number of relations
    //    initial_reln_inds = get_random_zero_one_list_and_increment_seed(
    //        total_reln_available, num_relns_requested);
    //    print(
    //        "\nWARNING Reading a maximum of {} relns. This is {}% of the available relations.\n".format(
    //            num_relns_requested,
    //            (float(num_relns_requested) / total_reln_available) * 100));
  }

  fmpz_mat_t sp_reln_mat;
  // set num rows
  fmpz_mat_init(sp_reln_mat, num_relns_requested, num_prs);
  slong nrows = fmpz_mat_ncols(sp_reln_mat);
  fmpz_mat_t lp_reln_mat;
  fmpz_mat_init(lp_reln_mat, nrows, nf.dg + 1);
  vector<QuarticElem> reln_gens = vector<QuarticElem>(nrows);
  //  elems_mat = Matrix(nf, nrows=sp_reln_mat.nrows(), ncols=1, sparse=True);
  //  print("Num of rows={} and cols={}".format(sp_reln_mat.nrows(), sp_reln_mat.ncols()));

  f.clear();
  f.seekg(0);
  vector<string> split_res{};
  vector<string> split_comma_res{};
  int is_not_valid0 = -1;
  int is_not_valid1 = -1;
  fmpq_t coeff0;
  fmpq_init(coeff0);
  fmpq_t coeff1;
  fmpq_init(coeff1);
  nf_elem_t elem_prod0;
  nf_elem_t elem_prod1;
  nf_elem_t elem;
  nf_elem_init(elem_prod0, nf.antic_nf);
  nf_elem_init(elem_prod1, nf.antic_nf);
  nf_elem_init(elem, nf.antic_nf);

  int idl_idx = 0;
  fmpz_t com;
  fmpz_init(com);
  fmpz_set_ui(com, idl_coeff_gcds[idl_idx]);
  nf_elem_t gen;
  nf_elem_init(gen, nf.antic_nf);
  nf_elem_set(gen, idl_poly_gens[idl_idx].elem, nf.antic_nf);
//  nf_elem_print_pretty(gen, nf.antic_nf, "x");
  nf_elem_t least_int;
  nf_elem_init(least_int, nf.antic_nf);
  nf_elem_set_ui(least_int, idl_least_ints[idl_idx], nf.antic_nf);

  //TODO
//  int loaded_reln_ind;
  std::pair <int, int> coeffs;
  set<pair<int, int>> seen;

  int32_t count_dups = 0;

  int reln_ind = 0;

  // 2nd pass over file
  set<int32_t> large_prs {};
  set<int32_t> small_prs {};
  for (i=0; i<max_num_lines; i++) {
    getline(f, line);
    if (f.eof()) {
      break;
    }
    if (StartsWith("#", line) or StartsWith(idl_poly_gen_str, line) or StartsWith(idl_coeff_gcd_str, line) or StartsWith(idl_fb_exp_vec_str, line)) {
      // Case: skip. This line does not trigger a change of ideal and is not a relation
      continue;
    } else if (StartsWith(idl_least_int_str, line)) {
      // Case: switch to a new ideal
      idl_idx += 1;
      if (idl_idx == idl_least_ints.size()) {
        // Take this to mean that we're at the end of the file
        break;
      }
      fmpz_set_ui(com, idl_coeff_gcds[idl_idx]);
      nf_elem_set(gen, idl_poly_gens[idl_idx].elem, nf.antic_nf);
      // the ideal we are sieving could be defined here
      continue;
    }

    //  set idl_poly_gen, coeff_gcd, idl_least_int
    split_res.clear();
    // get the coeffs a,b by splitting the string at ":"
    Split(string(":"), line, split_res);
    if (split_res.size() != 3) {
      string err_str = string("wrong number of ':' in line at index ") + to_string(i);
      throw runtime_error(err_str);
    }
    split_comma_res.clear();
    Split(string(","), split_res[0], split_comma_res);

    // remove relations that are a ZZ multiple of previous relation
    coeffs.first = stoi(split_comma_res[0]);
    coeffs.second = stoi(split_comma_res[1]);
    if (seen.count(coeffs)) {
      count_dups += 1;
      continue;
    } else {
      seen.insert(coeffs);
    }

    is_not_valid0 = fmpq_set_str(coeff0, split_comma_res[0].c_str(), 10);
    if (is_not_valid0) {
      string err_str = string("Failed getting fmpq from line at index ") + to_string(i) + "coeff0 " + split_comma_res[0];
      throw runtime_error(err_str);
    }
    is_not_valid1 = fmpq_set_str(coeff1, split_comma_res[1].c_str(), 10);
    if (is_not_valid1) {
      string err_str = string("Failed getting fmpq from line at index ") + to_string(i) + "coeff1 " + split_comma_res[1];
      throw runtime_error(err_str);
    }

    norm_pr_facs.clear();
    sieve_pr_facs_hex.clear();
    fmpz_t pr;
    fmpz_init(pr);
    Split(string(","), split_res[1], sieve_pr_facs_hex);
    for (auto hex_pr : sieve_pr_facs_hex) {
      norm_pr_facs.push_back(stoul(hex_pr, nullptr, 16));
    }
    large_prs.clear();
    small_prs.clear();
    fmpz_t prod_all_prs;
    fmpz_one(prod_all_prs);
    for (auto long_pr : norm_pr_facs) {
      fmpz_set_ui(pr, long_pr);
      fmpz_mul(prod_all_prs, prod_all_prs, pr);
      small_prs.insert(long_pr);
      if (long_pr > smoothness_bound) {
        large_prs.insert(long_pr);
      }
    }
    if (large_prs.size() > 1) {
      printf("More than one large prime divides the norm. Skippig this relation");
      continue;
    }
    fmpz_t elem_norm;
    fmpz_init(elem_norm);
    fmpz_mul_ui(elem_norm, com, idl_least_ints[idl_idx]);
    fmpz_mul(elem_norm, elem_norm, prod_all_prs);
    nf_elem_scalar_mul_fmpq(elem_prod0, least_int, coeff0, nf.antic_nf);
    nf_elem_scalar_mul_fmpq(elem_prod1, gen, coeff1, nf.antic_nf);
    nf_elem_set(elem, elem_prod0, nf.antic_nf);
    nf_elem_add(elem, elem, elem_prod1, nf.antic_nf);
    reln_gens.emplace_back(elem);

#ifndef NDEBUG
    fmpz_t elem_norm_expected;
    fmpz_init(elem_norm_expected);
    nf_elem_norm_abs(elem_norm_expected, elem, nf.antic_nf);
    assert(fmpz_equal(elem_norm, elem_norm_expected));
#endif

    fmpz_factor_t com_facs;
    fmpz_factor_init(com_facs);
    fmpz_factor_smooth(com_facs, com, lpb, false);
    for (int j=0; j<com_facs[0].num; j++) {
      small_prs.insert((com_facs[0].p)[j]);
    }
    fmpz_t least_int_fmpz;
//    fmpz_init(least_int);
    fmpz_init_set_ui(least_int_fmpz, idl_least_ints[idl_idx]);
    fmpz_factor_t idl_least_int_facs;
    fmpz_factor_init(idl_least_int_facs);
    fmpz_factor_smooth(idl_least_int_facs, least_int_fmpz, lpb, false);
    for (int j=0; j<idl_least_int_facs[0].num; j++) {
      small_prs.insert((idl_least_int_facs[0].p)[j]);
    }

    for (auto lp: large_prs) {
      // construct ideals for large primes
      nmod_poly_t mod_poly;
      nmod_poly_init(mod_poly, lp);
      _nmod_poly_set_length(mod_poly, 5);
      for (int k=0; k<nf.dg; k++) {
        cout << "assigning" << nf.poly_coeffs[k] << endl;
        nmod_poly_set_coeff_ui(mod_poly, nf.poly_coeffs[k], k);
      }
      nmod_poly_print(mod_poly);
//      poly_str_to_nmod_poly(mod_poly, poly_c_str);

//      nmod_poly_factor_t poly_facs;
//      nmod_poly_factor_init(poly_facs);
//      nmod_poly_factor(poly_facs, mod_poly);
//      nmod_poly_factor_print(poly_facs);

      // how to iterate over this ?

      // iterate
//      long res_class_dg =
//
//      // is there a way to convert fmpq poly to fmpz poly
//
//      fmpz_t lp_fmpz;
//      fmpz_init_set_ui(lp_fmpz, lp);
//      QuarticIdeal(, fmpq_poly_t g1_pol, fmpz_mat_t mat, bool is_prime, uint res_class_deg_in, uint ram_ind_in, QuarticField &nf)
    }
    reln_ind +=1;
  }
  cout << "num of lp" << large_prs.size() << endl;
  cout << reln_ind << endl;
  // compute the valuation of the elem
  cout << "number of dups is" << endl;
  cout << count_dups << endl;

  //delete elems in the array
//  for (auto _ : reln_gens) {
//    _.QuarticElemClear();
//  }
//  for (auto _ : idl_poly_gens) {
//    _.QuarticElemClear();
//  }
  f.close();
}
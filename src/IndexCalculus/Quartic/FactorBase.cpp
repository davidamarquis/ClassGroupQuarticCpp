#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

#include <antic/nf.h>
#include <antic/nf_elem.h>
#include <flint/arith.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include "FilterUtils.hpp"
#include "include/rapidjson/document.h"
#include "include/rapidjson/istreamwrapper.h"
#include "include/rapidjson/ostreamwrapper.h"
#include "include/rapidjson/stringbuffer.h"
#include "include/rapidjson/writer.h"
#include "src/IndexCalculus/Quartic/ObjStringConverters.hpp"
#include "src/Quartic/PredefinedIdeals.hpp"
#include "src/Quartic/QuarticField.hpp"
#include "src/Quartic/QuarticIdeal.hpp"
#include "src/Flint.hpp"
#include "FactorBase.hpp"
#include "src/IndexCalculus/Quartic/IoUtils.hpp"

#include <vector>
#include <inttypes.h>
#include <cstdlib>

using namespace std;
using namespace rapidjson;

FactorBase::FactorBase(vector<QuarticIdeal> prime_idls, vector<uint32_t> rat_primes, vector<uint8_t> num_primes_above, uint32_t fb_bound) : prime_idls(prime_idls), rat_primes(rat_primes), num_primes_above(num_primes_above), fb_bound(fb_bound) {}

inline const char* mat_idl_key = "mat";
inline const char* gen0_idl_key = "gen0";
inline const char* gen1_idl_key = "gen1";
inline const char* residue_class_degree_idl_key = "residue_class_degree";
inline const char* ramification_index_idl_key = "ramification_index";

inline const char* prime_idls_fb_key = "prime_idls";
inline const char* smoothness_bound_fb_key = "smoothness_bound";
inline const char* num_facs_above_prime_fb_key = "num_facs_above_prime";
inline const char* rat_primes_fb_key = "rat_primes";

void FactorBase::print() {
  if (prime_idls.size()==0) {
    printf("Factor base containing 0 prime ideals.");
  }
  printf("Factor Base:{Size %lu\nNF\n", prime_idls.size());
  prime_idls[0].nf.print();
  printf("Smoothness bound %" PRIu32 "}\n", fb_bound);
}

FactorBase FactorBaseFromJson(const char* fb_json_path, QuarticField &nf, uint32_t smoothness_bound) {
  Document fb_doc {};
  LoadDoc(fb_doc, fb_json_path);

  GenericArray prime_idls_strs = fb_doc[prime_idls_fb_key].GetArray();

  size_t num_elems_in_fb = static_cast<size_t>(prime_idls_strs.Size());
  string mat_str = "";
  vector<QuarticIdeal> prime_idls = {};
  prime_idls.reserve(num_elems_in_fb);

  fmpz_t g0;
  fmpq_poly_t g1;
  fmpz_mat_t mat;
  fmpz_init(g0);
  fmpq_poly_init(g1);
  fmpz_mat_init(mat, 4, 4);
  int res_class_deg;
  int ram_ind;
  fmpz_t nrm;
  fmpz_init(nrm);
  for (uint i=0;i<num_elems_in_fb;i++) {
    mat_str_to_fmpz_mat(mat, prime_idls_strs[i][mat_idl_key].GetString());
    fmpz_set_str(g0, prime_idls_strs[i][gen0_idl_key].GetString(), 10);
    assert(fmpz_is_prime(g0));
    // stop reading ideals when we hit the smoothness bound
    if (fmpz_cmp_ui(g0, smoothness_bound) > 0) {
      break;
    }
    poly_str_to_fmpq_poly(g1, prime_idls_strs[i][gen1_idl_key].GetString());
    res_class_deg = prime_idls_strs[i][residue_class_degree_idl_key].GetInt();
    ram_ind = prime_idls_strs[i][ramification_index_idl_key].GetInt();
    prime_idls.emplace_back(g0, g1, mat, true, res_class_deg, ram_ind, nf);
#ifndef NDEBUG
    fmpz_t idl_nrm;
    fmpz_init(idl_nrm);
    fmpz_t radical_of_nrm;
    fmpz_init(radical_of_nrm);

    if (!fmpz_poly_is_zero(g1)) {
      nf_elem_norm_abs(nrm, prime_idls.back().gen1, nf.antic_nf);
      fmpz_gcd(idl_nrm, nrm, g0);
      fmpz_is_perfect_power(radical_of_nrm, idl_nrm);
      bool is_nrm_prime = fmpz_is_prime(idl_nrm);
      if (not (fmpz_is_prime(radical_of_nrm) or is_nrm_prime)) {
        throw invalid_argument(string("prime ideal has non-prime radical. Index=") + to_string(i));
      }
    }
#endif
  }
  fmpz_clear(g0);
  fmpq_poly_clear(g1);
  fmpz_mat_clear(mat);

  vector<uint32_t> rat_primes {};
  vector<uint8_t> num_primes_above {};
  uint32_t fb_primes_bound = fb_doc[smoothness_bound_fb_key].GetInt();
  if (smoothness_bound > fb_primes_bound) {
    throw invalid_argument("smoothness_bound must be smaller than fb_primes_bound");
  }

  GenericArray rat_primes_strs = fb_doc[rat_primes_fb_key].GetArray();
  for (SizeType i = 0; i < rat_primes_strs.Size(); i++) {
    rat_primes.push_back(strtoull(rat_primes_strs[i].GetString(), nullptr, 10));
  }
  GenericArray num_primes_above_strs = fb_doc[num_facs_above_prime_fb_key].GetArray();
  for (SizeType i = 0; i < num_primes_above_strs.Size(); i++) {
    num_primes_above.push_back(strtoull(num_primes_above_strs[i].GetString(), nullptr, 10));
  }

  return FactorBase(prime_idls, rat_primes, num_primes_above, smoothness_bound);
}

void FactorBase::TrivialRelations(fmpz_mat_t triv_reln_mat) {
  uint32_t nrows = rat_primes.size();
  uint32_t ncols = prime_idls.size();
  if (fmpz_mat_nrows(triv_reln_mat) != nrows) {
    throw runtime_error("TrivialRelations: number of rows in matrix is not equal to the number of rational primes");
  }
  if (fmpz_mat_ncols(triv_reln_mat) != ncols) {
    throw runtime_error("TrivialRelations: number of cols in matrix is not equal to the number of prime ideals");
  }
  uint32_t col_ind = 0;
  fmpz_t entry_value;
  fmpz_init(entry_value);

  fmpz_mat_init(triv_reln_mat, nrows, ncols);
  for(uint32_t rr=0;rr<nrows;rr++) {
    for(uint32_t cc=col_ind;cc<(col_ind + num_primes_above[rr]);cc++) {
      fmpz_set_ui(entry_value, prime_idls[cc].ram_ind);
      fmpz_set(fmpz_mat_entry(triv_reln_mat, rr, cc), entry_value);
    }
    col_ind += num_primes_above[rr];
  }
  fmpz_clear(entry_value);
}


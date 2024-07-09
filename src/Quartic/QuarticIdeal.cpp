#include "QuarticIdeal.hpp"

#include <cassert>
#include <algorithm>
#include <flint/fmpz.h>
#include <src/Flint.hpp>

void QuarticIdeal::InitValueMemberVariables() {
  fmpz_init(gen0);
  nf_elem_init(gen1, nf.antic_nf);
  fmpz_mat_init(basis_mat, 4, 4);
  nf_elem_init(b0, nf.antic_nf);
  nf_elem_init(b1, nf.antic_nf);
  nf_elem_init(b2, nf.antic_nf);
  nf_elem_init(b3, nf.antic_nf);
  nf_elem_init(tau0, nf.antic_nf);
  nf_elem_init(uniformizer, nf.antic_nf);
}

void QuarticIdeal::Assign(const QuarticIdeal &idl) {
  is_prime = idl.is_prime;

  res_class_deg = idl.res_class_deg;
  ram_ind = idl.ram_ind;
  fmpz_set(gen0, idl.gen0);
  nf_elem_set(gen1, idl.gen1, nf.antic_nf);
  fmpz_mat_set(basis_mat, idl.basis_mat);
  nf_elem_set(b0, idl.b0, nf.antic_nf);
  nf_elem_set(b1, idl.b1, nf.antic_nf);
  nf_elem_set(b2, idl.b2, nf.antic_nf);
  nf_elem_set(b3, idl.b3, nf.antic_nf);
  nf_elem_set(tau0, idl.tau0, nf.antic_nf);
  nf_elem_set(uniformizer, idl.uniformizer, nf.antic_nf);
}

QuarticIdeal::QuarticIdeal(QuarticField &nf) : nf(nf) {
  InitValueMemberVariables();
}

QuarticIdeal::QuarticIdeal(const QuarticIdeal&& idl) : nf(idl.nf) {
  InitValueMemberVariables();
  Assign(idl);
}

QuarticIdeal::QuarticIdeal(const QuarticIdeal& idl) : nf(idl.nf) {
  InitValueMemberVariables();
  Assign(idl);
}

QuarticIdeal::QuarticIdeal(fmpz_t g0, fmpq_poly_t g1_pol, fmpz_mat_t mat, bool is_prime, uint res_class_deg_in, uint ram_ind_in, QuarticField &nf) : is_prime(is_prime), res_class_deg(res_class_deg_in), ram_ind(ram_ind_in), nf(nf) {
  if (res_class_deg_in <= 0) {
    throw std::invalid_argument("Residue class degree must be positive");
  }
  if (ram_ind_in <= 0) {
    throw std::invalid_argument("Ramification index must be positive");
  }
  if (fmpz_mat_nrows(mat) != 4 || fmpz_mat_ncols(mat) != 4)  {
    throw std::invalid_argument("Mat is not a 4x4 matrix");
  }

  if (!(is_prime && fmpz_is_probabprime(g0))) {
    throw std::invalid_argument("Prime ideals about a rational prime pr must have g0=pr");
  }
  InitValueMemberVariables();

  fmpz_set(gen0, g0);
  nf_elem_set_fmpq_poly(gen1, g1_pol, nf.antic_nf);
  fmpz_mat_set(basis_mat, mat);
  BasisRowMulBasisElems(b0, 0);
  BasisRowMulBasisElems(b1, 1);
  BasisRowMulBasisElems(b2, 2);
  BasisRowMulBasisElems(b3, 3);

  // elements used for valuations are initialized to 0 and will be lazily initialized later
  nf_elem_set_si(tau0, 1, nf.antic_nf);
  nf_elem_set_si(uniformizer, 1, nf.antic_nf);
}

QuarticIdeal& QuarticIdeal::operator=(const QuarticIdeal &other) {
  if (this != &other) {
    Assign(other);
  }
  return *this;
}

QuarticIdeal& QuarticIdeal::operator=(const QuarticIdeal &&other) {
  if (this != &other) {
    Assign(other);
  }
  return *this;
}

bool QuarticIdeal::operator==(const QuarticIdeal& other) {
  /*
   * this operator tests if the two ideals are represented in the same way
   */
  bool eq_gen0 = fmpz_equal(gen0, other.gen0);
  bool eq_gen1 = nf_elem_equal(gen1, other.gen1, nf.antic_nf);
  bool eq_mat = fmpz_mat_equal(basis_mat, other.basis_mat);
  return (nf == other.nf) & eq_gen0 & eq_gen1 & eq_mat;
}

QuarticIdeal::~QuarticIdeal() {
//  printf("destroy idl w/ gen0, gen1 ");
//  fmpz_print(gen0);
//  printf(" ");
//  nf_elem_print_pretty(gen1, nf.antic_nf, "x");
//  printf("\n");

  fmpz_clear(gen0);
  nf_elem_clear(gen1, nf.antic_nf);
  fmpz_mat_clear(basis_mat);
  nf_elem_clear(b0, nf.antic_nf);
  nf_elem_clear(b1, nf.antic_nf);
  nf_elem_clear(b2, nf.antic_nf);
  nf_elem_clear(b3, nf.antic_nf);
  nf_elem_clear(tau0, nf.antic_nf);
  nf_elem_clear(uniformizer, nf.antic_nf);
}

void QuarticIdeal::GetBasisElem(int i, nf_elem_t &bas_elem) {
  if (i==0) {
    nf_elem_set(bas_elem, b0, nf.antic_nf);
  } else if (i==1) {
    nf_elem_set(bas_elem, b1, nf.antic_nf);
  } else if (i==2) {
    nf_elem_set(bas_elem, b2, nf.antic_nf);
  } else if (i==3) {
    nf_elem_set(bas_elem, b3, nf.antic_nf);
  } else {
    assert(false);
  }
}

void QuarticIdeal::Fprint(FILE *fp) {
  fmpq_poly_t pol;
  fmpq_poly_init(pol);
  nf_elem_get_fmpq_poly(pol, gen1, nf.antic_nf);

  fprintf(fp, "(");
  fmpz_fprint(fp, gen0);
  fprintf(fp, ",");
  fmpq_poly_fprint_pretty(fp, pol,"x");
  fprintf(fp, ")\n");

  fmpq_poly_clear(pol);
}

void QuarticIdeal::Print() { Fprint(stdout);
}

void QuarticIdeal::RandomElem(nf_elem_t elem, flint_rand_t state) {
  if (fmpz_mat_is_zero(basis_mat)) {
    RandomElemWithGens(elem, state);
    return;
  }

  RandomElemWithZBasis(elem, state);
}

void QuarticIdeal::RandomElemWithZBasis(nf_elem_t elem, flint_rand_t state) {
  fmpz_t c0;
  fmpz_t c1;
  fmpz_t c2;
  fmpz_t c3;
  fmpz_init(c0);
  fmpz_init(c1);
  fmpz_init(c2);
  fmpz_init(c3);
  fmpz_randm(c0, state, gen0);
  fmpz_randm(c1, state, gen0);
  fmpz_randm(c2, state, gen0);
  fmpz_randm(c3, state, gen0);

  nf_elem_t prd0;
  nf_elem_init(prd0, nf.antic_nf);
  nf_elem_scalar_mul_fmpz(prd0, b0, c0, nf.antic_nf);
  nf_elem_t prd1;
  nf_elem_init(prd1, nf.antic_nf);
  nf_elem_scalar_mul_fmpz(prd1, b1, c1, nf.antic_nf);
  nf_elem_t prd2;
  nf_elem_init(prd2, nf.antic_nf);
  nf_elem_scalar_mul_fmpz(prd2, b2, c2, nf.antic_nf);
  nf_elem_t prd3;
  nf_elem_init(prd3, nf.antic_nf);
  nf_elem_scalar_mul_fmpz(prd3, b3, c3, nf.antic_nf);

  nf_elem_add(elem, elem, prd0, nf.antic_nf);
  nf_elem_add(elem, elem, prd1, nf.antic_nf);
  nf_elem_add(elem, elem, prd2, nf.antic_nf);
  nf_elem_add(elem, elem, prd3, nf.antic_nf);
  nf_elem_reduce(elem, nf.antic_nf);

  fmpz_clear(c0);
  fmpz_clear(c1);
  fmpz_clear(c2);
  fmpz_clear(c3);
  nf_elem_clear(prd0, nf.antic_nf);
  nf_elem_clear(prd1, nf.antic_nf);
  nf_elem_clear(prd2, nf.antic_nf);
  nf_elem_clear(prd3, nf.antic_nf);
}

void QuarticIdeal::RandomElemWithGens(nf_elem_struct* elem,
                                      flint_rand_s* state) {
  /*
   * A fallback function in case we can't call RandomElemWithZBasis
   */
  nf_elem_t r0;
  nf_elem_init(r0, nf.antic_nf);
  nf_elem_randtest(r0, state, 10, nf.antic_nf);

  nf_elem_t r1;
  nf_elem_init(r1, nf.antic_nf);
  nf_elem_randtest(r1, state, 10, nf.antic_nf);

  nf_elem_scalar_mul_fmpz(r0, r0, gen0, nf.antic_nf);
  nf_elem_mul(r1, r1, gen1, nf.antic_nf);

  nf_elem_add(elem, r0, r1, nf.antic_nf);
}

int QuarticIdeal::Valuation(nf_elem_t elem) {
  /*
   * Returns the valuation of an elem of Ok. Result is undefined if elem is not in Ok
   */
#ifndef NDEBUG
  fmpz_mat_t elem_coords;
  fmpz_mat_init(elem_coords,1,4);
  if(!nf.ElemToOkCoordinates(elem_coords, elem)) {
    throw std::invalid_argument("Elem must be in Ok");
  }
#endif
  // this function will silently fail if elem does not have the same nf as self
  if (nf_elem_is_one(tau0, nf.antic_nf)) {
    throw std::invalid_argument("InitTau0 must be called before this function");
  }
  if (!is_prime) {
    throw std::invalid_argument("Valuation only defined for prime ideals");
  }
  nf_elem_t yy; // an elem of Ok
  nf_elem_init(yy, nf.antic_nf);
  int nrm_valn = nf_elem_norm_valn(gen0, elem, nf.antic_nf);
  int max_valn;
  max_valn = max(1, (int)ceil(nrm_valn / res_class_deg));
  int valn = 0;

  nf_elem_t y_prime; // an elem of Ok
  nf_elem_init(y_prime, nf.antic_nf);

  nf_elem_set(yy, elem, nf.antic_nf);

  fmpz_mat_t coords;
  fmpq_mat_t coords_q;
  fmpz_mat_init(coords, 1, 4);
  fmpq_mat_init(coords_q, 1, 4);
  slong pr = fmpz_get_si(gen0);
  fmpz_t dnm_std_coords;
  fmpz_init(dnm_std_coords);
  fmpz_t dnm_Ok_coords;
  fmpz_init(dnm_Ok_coords);

  nmod_mat_t coords_Fp;
  nmod_mat_init(coords_Fp, 1, 4, pr);
  bool ret;
  for (int k=0;k<=max_valn;k++) {
    nf_elem_mul(y_prime, yy, tau0, nf.antic_nf);
    ret = nf.ElemToOkCoordinates(coords, y_prime);
    if (!ret) {
      throw runtime_error("Encountered elem that is not in Ok");
    }
    fmpz_mat_get_nmod_mat(coords_Fp, coords);
    if (!nmod_mat_is_zero_row(coords_Fp, 0)) {
      break;
    }
    nf_elem_scalar_div_si(yy, y_prime, pr, nf.antic_nf);
    valn += 1;
  }
  nf_elem_clear(yy, nf.antic_nf);
  nf_elem_clear(y_prime, nf.antic_nf);
  fmpz_mat_clear(coords);
  fmpq_mat_clear(coords_q);
  nmod_mat_clear(coords_Fp);
  fmpz_clear(dnm_std_coords);
  fmpz_clear(dnm_Ok_coords);
  return valn;
}

void QuarticIdeal::InitTau0(fmpq_mat_t uniformizer_mul_mat, flint_rand_t state) {
  /*
   * Performs lazy initialization of the member variable tau0
   * tau0 is a field element so that tau0*uniformizer = 0 (mod gen0). From Belabas prop 5.9
   */
  // proceed only if this function has not already run
  if (!nf_elem_is_one(tau0, nf.antic_nf)) {
    return;
  }

  if (!is_prime) {
    throw std::invalid_argument("Tau0 is only defined for prime ideals");
  }

  if (!nf_elem_is_zero(uniformizer, nf.antic_nf)) {
    InitUniformizer(state);
  }

  fmpz_mat_t bas_coords_z;
  fmpz_mat_init(bas_coords_z,4,4);

  nf_elem_t elem;
  nf_elem_init(elem, nf.antic_nf);
  fmpz_t dnms[4];
  for (int i=0;i<4;i++) {
    fmpz_init(dnms[i]);
  }

  // construct coordinate matrix of multiplication of uniformizer and the integral basis
  for (int i=0;i<4;i++) {
    nf.GetBasisElem(i, elem);
    nf_elem_mul(elem, elem, uniformizer, nf.antic_nf);
    nf_elem_get_fmpz_mat_row(bas_coords_z, i, dnms[i], elem, nf.antic_nf);
  }
  fmpq_mat_set_fmpz_mat(uniformizer_mul_mat, bas_coords_z);
  fmpq_mat_divide_rows(uniformizer_mul_mat, dnms);
  fmpq_mat_mul(uniformizer_mul_mat, uniformizer_mul_mat, nf.int_basis_mat_inv);

  fmpz_mat_t uniformizer_mul_mat_num;
  fmpz_mat_init(uniformizer_mul_mat_num,4,4);
  fmpz_t uniformizer_mul_mat_dnm;
  fmpz_init(uniformizer_mul_mat_dnm);
  fmpq_mat_get_fmpz_mat_matwise(uniformizer_mul_mat_num, uniformizer_mul_mat_dnm, uniformizer_mul_mat);

  // init the matrix (mod p)
  nmod_mat_t uniformizer_mul_mod_mat;
  nmod_mat_init(uniformizer_mul_mod_mat,4,4, fmpz_get_si(gen0));
  fmpz_mat_get_nmod_mat(uniformizer_mul_mod_mat, uniformizer_mul_mat_num);

  // get the kernel
  nmod_mat_transpose(uniformizer_mul_mod_mat, uniformizer_mul_mod_mat);
  nmod_mat_t kern;
  nmod_mat_init(kern,4,4,fmpz_get_si(gen0));
  ulong rank = nmod_mat_nullspace(kern, uniformizer_mul_mod_mat);
  if (rank == 0) {
    throw std::length_error("Nullspace should contain at least one vector");
  }
  nmod_mat_transpose(kern, kern);

  fmpz_mat_t kern_z;
  fmpz_mat_init(kern_z,4,4);
  fmpz_mat_set_nmod_mat(kern_z, kern);
  fmpz_t f_one;
  fmpz_init(f_one);
  fmpz_one(f_one);

  nf.OkCoordinateMatRowToElem(tau0, kern_z,
                              0); // assign tau0 to elem corresponding to the first vector in the kernel

  fmpz_clear(f_one);
  fmpz_mat_clear(kern_z);
  nmod_mat_clear(kern);
  nmod_mat_clear(uniformizer_mul_mod_mat);
  fmpz_mat_clear(uniformizer_mul_mat_num);
  fmpz_clear(uniformizer_mul_mat_dnm);
  for (int i=0;i<4;i++) {
    fmpz_clear(dnms[i]);
  }
  fmpz_mat_clear(bas_coords_z);
  nf_elem_clear(elem, nf.antic_nf);
}

void QuarticIdeal::InitUniformizer(flint_rand_t state) {
  /*
   * Performs lazy initialization of the member variable uniformizer
   * this randomized algorithm is Belabas 6.1.1
  */
  // proceed only if this function has not already run
  if (!nf_elem_is_one(uniformizer, nf.antic_nf)) {
    return;
  }

  if (!is_prime) {
    throw std::invalid_argument("Uniformizer is only defined for prime ideals");
  }
  nf_elem_t rand;
  nf_elem_init(rand, nf.antic_nf);

  int max_tries = 1000; // the expected number of tries is <= 1/(1-1/p)^4 for an ideal above the prime p. Worst case is 2 completely splits where the expected number of tries is 16
  slong ss = 0;
  int i = 0;
  for(i=0;i<max_tries; i++) {
    RandomElem(rand, state);
    ss = nf_elem_norm_valn(gen0, rand, nf.antic_nf);
    assert(ss >= 1);
    if (ss < res_class_deg + 1) {
      break;
    }
  }
  assert(ss < res_class_deg + 1);
  if (i == max_tries+1) {
    throw std::invalid_argument("Failed to find uniformizer");
  }
  nf_elem_set(uniformizer, rand, nf.antic_nf);
  nf_elem_clear(rand, nf.antic_nf);
}

void QuarticIdeal::BasisRowMulBasisElems(nf_elem_t new_elem, const slong r) {
  nf_elem_t tmp_elem;
  nf_elem_t bb;
  nf_elem_init(tmp_elem, nf.antic_nf);
  nf_elem_init(bb, nf.antic_nf);
  for (slong cc=0; cc<fmpz_mat_ncols(basis_mat); cc++) {
    nf_elem_set_fmpz(tmp_elem, fmpz_mat_entry(basis_mat, r, cc), nf.antic_nf);
    nf.GetBasisElem(cc, bb);
    nf_elem_mul(tmp_elem, tmp_elem, bb, nf.antic_nf);
    nf_elem_add(new_elem, new_elem, tmp_elem, nf.antic_nf);
  }
  nf_elem_clear(tmp_elem, nf.antic_nf);
  nf_elem_clear(bb, nf.antic_nf);
}

void QuarticIdealFromCParams(QuarticIdeal& idl, QuarticField& nf, const ulong gen0_u, const char* gen1_str, const char* mat_str, uint res_class_deg, uint ram_ind) {
  fmpz_t g0;
  fmpz_init_set_ui(g0, gen0_u);
  fmpq_poly_t g1;
  fmpq_poly_init(g1);
  poly_str_to_fmpq_poly(g1, gen1_str);

  fmpz_mat_t mat;
  fmpz_mat_init(mat,4,4);
  if (mat_str != nullptr) {
    mat_str_to_fmpz_mat(mat, mat_str);
  } else {
    mat_str = "1, 0, 0, 0\n0, 1, 0, 0\n0, 0, 1, 0\n0, 0, 0, 1";
  }

  idl = QuarticIdeal(g0, g1, mat, true, res_class_deg, ram_ind, nf);
  fmpz_clear(g0);
  fmpq_poly_clear(g1);
  fmpz_mat_clear(mat);
}

#ifndef THESIS_QUARTICIDEAL_HPP
#define THESIS_QUARTICIDEAL_HPP

#include <src/Quartic/QuarticField.hpp>
class QuarticField;

class QuarticIdeal {
  public:
  bool is_prime;
  int res_class_deg;
  int ram_ind;
  fmpz_t gen0; //gen0 is the smallest integer in the ideal
  nf_elem_t gen1; //gen1 is chosen so that (gen0, gen1) is a 2-element representation
  fmpz_mat_t basis_mat;
  nf_elem_t b0;
  nf_elem_t b1;
  nf_elem_t b2;
  nf_elem_t b3;
  QuarticField &nf;

  nf_elem_t uniformizer;
  nf_elem_t tau0;

  void InitValueMemberVariables();

  void Assign(const QuarticIdeal &idl);

  QuarticIdeal() = delete;

  QuarticIdeal(QuarticField &nf);

  QuarticIdeal(const QuarticIdeal&& idl);

  QuarticIdeal(const QuarticIdeal& idl);

  QuarticIdeal(fmpz_t g0, fmpq_poly_t g1_pol, fmpz_mat_t mat, bool is_prime, uint res_class_deg_in, uint ram_ind_in, QuarticField &nf);

  QuarticIdeal& operator=(const QuarticIdeal &other);

  QuarticIdeal& operator=(const QuarticIdeal &&other);

  bool operator==(const QuarticIdeal& other);

  ~QuarticIdeal();

  void GetBasisElem(int i, nf_elem_t &bas_elem);

  void Fprint(FILE *fp);

  void Print();

  void RandomElem(nf_elem_t elem, flint_rand_t state);
  void RandomElemWithGens(nf_elem_t elem, flint_rand_t state);
  void RandomElemWithZBasis(nf_elem_t elem, flint_rand_t state);

  int Valuation(nf_elem_t elem);
  void InitUniformizer(flint_rand_t state);
  void InitTau0(fmpq_mat_t uniformizer_mul_mat, flint_rand_t state);

  private:
  void BasisRowMulBasisElems(nf_elem_t new_elem, const slong r);
};

void QuarticIdealFromCParams(QuarticIdeal &idl, QuarticField &nf, const ulong gen0_u, const char* gen1_str, const char* mat_str, uint res_class_deg, uint ram_ind);

#endif  // THESIS_QUARTICIDEAL_HPP
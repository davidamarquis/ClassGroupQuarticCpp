#ifndef THESIS_FACTORBASE_HPP
#define THESIS_FACTORBASE_HPP

#include <cstdint>
#include <vector>
#include <src/Quartic/QuarticIdeal.hpp>
class QuarticField;

class FactorBase {
  public:
  explicit FactorBase(vector<QuarticIdeal> prime_idls, vector<uint32_t> rat_primes, vector<uint8_t> num_primes_above, uint32_t smoothness_bound);
  FactorBase() = default;
  FactorBase(const FactorBase &) = default;
  FactorBase(FactorBase &&) = default;
  FactorBase& operator=(const FactorBase &) = default;
  FactorBase& operator=(FactorBase &&) = default;
  ~FactorBase() = default;

  vector<QuarticIdeal> prime_idls;
  vector<uint32_t> rat_primes;
  vector<uint8_t> num_primes_above;
  uint32_t fb_bound;

  void print();

  void TrivialRelations(fmpz_mat_t triv_reln_mat);
};

FactorBase FactorBaseFromJson(const char* fb_json_path, QuarticField &nf, uint32_t smoothness_bound);

#endif  // THESIS_FACTORBASE_HPP
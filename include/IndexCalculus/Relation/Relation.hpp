#ifndef RELATION_H
#define RELATION_H

#include <vector>
#include "ANTL/Interface/Multiplicative.hpp"
#include <ANTL/IndexCalculus/FactorBase/FactorBase.hpp>
#include <ANTL/common.hpp>

using namespace ANTL;

namespace ANTL {

  class Relation
  {
  protected:
  	static long sizeFB;         // factor base size (length of relations)
  	IMultiplicative gamma;  // the principal ideal generator

  	std::vector<long> vec_idx;
  	std::vector<long> vec_exp;

  public:

  	Relation() = default;

  	static void set_sizeFB(long newsize);
  	static long get_sizeFB();

  	std::vector<long> get_vec_idx() const {return vec_idx;}
  	std::vector<long> get_vec_exp() const {return vec_exp;}

  	long get_idx(long i) const {}
  	long get_exp(long i) const {}

  	void set_exp ( long i , long exp_in ) {vec_exp[i] = exp_in;}

  	long get_size() const {return (long)vec_idx.size();}

  	/* Adds a new factor to the list */
  	void add_element (long idx_in , long exp_in) {}

  	/* Takes out the factor of index idx_in from the
  	 * factor list.
  	 */
  	void del_element (long idx_in) {}

  	void assign_zero() {}

    void assign(const Relation &rel) {}
  	virtual void assign(const std::vector<long> & indices , const std::vector<long> & exponents, const IMultiplicative &q) {}
  	virtual void set_gamma(const IMultiplicative &q) { gamma=q; }

  	/* Checks out if the power product equals gamma */
  	int check(const FactorBase &fac_base) {return 0;}
  };

} // ANTL

#endif // guard

#ifndef QUADCLASSGROUPINDCALC_H
#define QUADCLASSGROUPINDCALC_H
#include <string>
#include <vector>
#include <memory>
#include <NTL/ZZ.h>
#include <ANTL/IndexCalculus/IndCalc/IndCalc.hpp>
#include <ANTL/Interface/OrderInvariants.hpp>
#include <ANTL/IndexCalculus/RelationGenerator/QuadRelationGenerator.hpp>
#include <ANTL/IndexCalculus/Relation/QuadRelation.hpp>
#include <ANTL/IndexCalculus/FactorBase/QuadFactorBase.hpp>

template <class T, class R>
class QuadIndCalc:public IndCalc<T,R> {
  /*TODO this class is an example to show which functions a class inheriting from IndCalc needs to implement */

public:
  // initialization
  static std::unique_ptr<QuadIndCalc<T,R>> create(IOrder<T,R> const &order, std::map<std::string, std::string> const &params) {
    // Because of the inheritance structure of IndCalc, we need to do a two step initialization.
    // This initializes an IndCalc with shared pointers for FactorBase and RelationGenerator as member variables
    unique_ptr<QuadFactorBase> fac_base{ new QuadFactorBase(order, params) };
    unique_ptr<QuadRelationGenerator> reln_generator{ new QuadRelationGenerator(order, params, fac_base.get())};
    unique_ptr<QuadIndCalc<T,R>> ind_calc {new QuadIndCalc<T,R>(std::move(fac_base), std::move(reln_generator))};
    ind_calc->setup_mat();
    return ind_calc;
  }

  friend void FactorBase::push_to_fb(IMultiplicative &fb_elem); // ind calc can add elems to the factor base

  //TODO fill in these stubs after we have implemented an index calculus algorithm for this class
  virtual NTL::ZZ class_number() {return NTL::ZZ(0);};
  virtual std::vector<NTL::ZZ> class_group() {std::vector<NTL::ZZ> cg = {NTL::ZZ(4), NTL::ZZ(5)}; return cg;};
  virtual std::vector<T> unit_group() {std::vector<T> ug = {T()}; return ug;};
  virtual R regulator() {R reg = {R()}; return reg;};

  // implement pure virtual functions from base class
  RelationGenerator* const get_relation_generator() override {return this->relation_generator.get();};
  FactorBase* const get_factor_base() override {return this->factor_base.get();};

  virtual ~QuadIndCalc() = default;
protected:
  // initialization
  using IndCalc<T,R>::IndCalc; // inherit the constructors
  QuadIndCalc<T,R>(std::unique_ptr<QuadFactorBase> factor_base, std::unique_ptr<QuadRelationGenerator> relation_generator) :
    factor_base(std::move(factor_base)),
    relation_generator(std::move(relation_generator)) {}

  QuadIndCalc<T,R>() = default;
  QuadIndCalc(const QuadIndCalc&) = delete;
  QuadIndCalc& operator=(const QuadIndCalc&) = delete;

  // these pure virtual function implementations are defined outside this class declaration
  void compute_fac_base() override;
  void compute_relations() override;
  void compute_mat() override;
private:
  std::unique_ptr<QuadFactorBase> factor_base;
  std::unique_ptr<QuadRelationGenerator> relation_generator;
};

// For now these 3 functions stubs. The code in them is just to suggest what they need to do

template <class T, class R>
void QuadIndCalc<T,R>::compute_fac_base() {
  IMultiplicative fb_elem = IMultiplicative();

  for(int i=0; i < get_factor_base()->get_size(); i++) {
    get_factor_base()->push_to_fb(fb_elem);
  }
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_relations() {
  bool result;
  /* compute_relations does sieving or random exponents to construct a relation */

  Relation quad_relation = Relation();
  auto reln_gen = get_relation_generator();
  for(long i=0; i < reln_gen->get_max_num_tests(); i++) {
    result = reln_gen->get_relation(quad_relation, i);
    if (result) {
      IndCalc<T,R>::relations.push_back(quad_relation);
    }
  }
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_mat() {
  //cf ANTL/include/ANTL/quadratic/index_calculus/qo_relation_matrix.hpp
  IndCalc<T,R>::rels_mat.SetDims(IndCalc<T,R>::relations.size(), get_factor_base()->get_size());

  // set matrix using functions in QuadRelation
};

#include "src/IndexCalculus/IndCalc/QuadIndCalc_impl.hpp"

#endif //QUADCLASSGROUPINDCALC_H

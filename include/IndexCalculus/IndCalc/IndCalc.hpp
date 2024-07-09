#ifndef CLASSGROUPINDCALC_H
#define CLASSGROUPINDCALC_H

#include <vector>
#include <string>
#include <map>
#include <memory>
#include <ANTL/common.hpp>
#include <ANTL/Interface/OrderInvariants.hpp>
#include <ANTL/IndexCalculus/FactorBase/FactorBase.hpp>
#include <ANTL/IndexCalculus/RelationGenerator/RelationGenerator.hpp>

using namespace ANTL;

template <class T, class R> // type of unit, type of regulator
class IndCalc : IOrder<T,R> {
public:
  std::vector<Relation> relations;
  NTL::Mat<ZZ> rels_mat;

  // subclasses must implement all four invariants class_number, class group, unit group, regulator
  virtual NTL::ZZ class_number() = 0;
  virtual std::vector<NTL::ZZ> class_group() = 0;
  virtual std::vector<T> unit_group() = 0;
  virtual R regulator() = 0;

  // subclasses need to choose the member variables that are returned by these getters
  virtual FactorBase* const get_factor_base() = 0;
  virtual RelationGenerator* const get_relation_generator() = 0;

  virtual ~IndCalc() = default;
protected:
  // derived class calls setup_mat to start the index calculus computation
  void setup_mat();
  // these functions are called by setup_mat and must be implemented by the subclasses
  virtual void compute_fac_base() = 0;
  virtual void compute_relations() = 0;
  virtual void compute_mat() = 0;

  IndCalc<T,R>() = default;
private:
  void setup_fac_base();
  void setup_relations();
};

template <class T, class R> void IndCalc<T,R>::setup_fac_base() {
  compute_fac_base();
  //TODO
//  if (factor_base.size() == 0) {
    // raise error
//  }
}

template <class T, class R> void IndCalc<T,R>::setup_relations() {
  setup_fac_base();
  compute_relations();
  //TODO
//  if (relations.size() == 0) {
    // raise error
//  }
}

template <class T, class R> void IndCalc<T,R>::setup_mat() {
  setup_relations();
  compute_mat();
  if (rels_mat.NumRows() == 0 or rels_mat.NumCols() == 0) {
    // raise error
  }
}
#endif //CLASSGROUPINDCALC_H

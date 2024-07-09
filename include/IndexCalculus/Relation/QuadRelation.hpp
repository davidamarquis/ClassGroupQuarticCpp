#ifndef RELATION_H
#define RELATION_H

#include "Relation.hpp"
#include "QuadFactorBase.hpp"

namespace ANTL {

  template<class T> // type of gamma
  class QuadRelation:public Relation
  {
  public:
  	QuadRelation<T>();

  	// get the minimum of a relation
//    QuadraticNumber<T> get_minimum() const;

//    int check(const QuadFactorBase<T>& fac_base);
  };
} // ANTL

#include "src/IndexCalculus/Relation/QuadRelation_impl.hpp"

#endif // guard

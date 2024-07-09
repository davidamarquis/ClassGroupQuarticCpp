#ifndef QUADFACTORBASE_H
#define QUADFACTORBASE_H

#include <ANTL/IndexCalculus/FactorBase/FactorBase.hpp>
#include <ANTL/Interface/OrderInvariants.hpp>

namespace ANTL {

  template < class T > class QuadraticOrder; // an order that inherits from IOrder

  class QuadFactorBase : public FactorBase {
  public:
    using FactorBase::FactorBase; // inherit the construtors
  };
}

#include "src/IndexCalculus/FactorBase/QuadFactorBase_impl.hpp"

#endif //QUADFACTORBASE_H


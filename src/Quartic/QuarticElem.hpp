#ifndef THESIS_QUARTICELEM_HPP
#define THESIS_QUARTICELEM_HPP

#include <antic/nf.h>
#include <antic/nf_elem.h>
#include "src/IndexCalculus/Quartic/QuarticCommon.hpp"

class QuarticElem {
  public:
  nf_elem_t elem;

  // give a copy constructor for this class
  QuarticElem(const QuarticElem& other) {
//    nf_elem_clear(elem, antic_nf);
    nf_elem_init(elem, antic_nf);
    nf_elem_set(elem, other.elem, antic_nf);
  }

  //q: does this class satisfy the rule of 5 ?

  // add functions to this class to make it satisfy the rule of 5
  QuarticElem& operator=(const QuarticElem& other) {
//    nf_elem_clear(elem, antic_nf);
//    nf_elem_init(elem, antic_nf);
    nf_elem_set(elem, other.elem, antic_nf);
    return *this;
  }
  QuarticElem(QuarticElem&& other) {
    nf_elem_init(elem, antic_nf);
    nf_elem_set(elem, other.elem, antic_nf);
  }
  QuarticElem& operator=(QuarticElem&& other) {
    nf_elem_set(elem, other.elem, antic_nf);
    return *this;
  }

  QuarticElem(nf_elem_t Elem) {
    nf_elem_init(elem, antic_nf);
    nf_elem_set(elem, Elem, antic_nf);
  };
  QuarticElem() {
    nf_elem_init(elem, antic_nf);
    nf_elem_one(elem, antic_nf);
  };
  ~QuarticElem() {
    nf_elem_clear(elem, antic_nf);
    //This is for container compatibility. Use QuarticElemClear to delete objects in this class
  };

//  void QuarticElemClear() {
//    nf_elem_clear(elem, antic_nf);
//  }
};

#endif  // THESIS_QUARTICELEM_HPP
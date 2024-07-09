#ifndef INDCALC_TEST
#define INDCALC_TEST

#include <string>
#include <map>
#include <iostream>
#include <NTL/ZZ.h>
#include "../catch.hpp"
#include "ANTL/IndexCalculus/IndCalc/QuadIndCalc.hpp"
#include <ANTL/Interface/OrderInvariants.hpp>
#include "ANTL/Constants.hpp"
#include "ANTL/Quadratic/QuadraticOrder.hpp"

using namespace Constants;
using namespace NTL;
using namespace ANTL;

std::map<std::string, std::string> get_params(std::string size_fb_str,
                                              std::string bound_fb_str,
                                              std::string max_num_tests_str) {
  std::map<std::string, std::string> params {{num_relations, "0"}, {size_fb, size_fb_str}, {bound_fb, bound_fb_str},
                                             {max_num_tests, max_num_tests_str}};
  return params;
}

TEST_CASE("IndCalc: Create initializes IndCalc with correct parameters", "[IndCalc]") {
  QuadraticOrder<ZZ> order = QuadraticOrder<ZZ>(ZZ(13));
  auto expected_size_fb1 = 0;
  auto expected_bound_fb1 = 0;
  auto expected_max_num_tests1 = 0;
  auto ind_calc1 = QuadIndCalc<ZZ, RR>::create(order, get_params("0", "0", "0"));
  REQUIRE(ind_calc1->get_relation_generator()->get_size_fb() == expected_size_fb1);
  REQUIRE(ind_calc1->get_factor_base()->get_size_fb() == expected_size_fb1);
  REQUIRE(ind_calc1->get_factor_base()->get_bound() == expected_bound_fb1);
  REQUIRE(ind_calc1->get_relation_generator()->get_max_num_tests() == expected_max_num_tests1);

  auto expected_size_fb2 = 1;
  auto expected_bound_fb2 = 2;
  auto expected_max_num_tests2 = 3;
  auto ind_calc2 = QuadIndCalc<ZZ, RR>::create(order, get_params("1", "2", "3"));
  REQUIRE(ind_calc2->get_relation_generator()->get_size_fb() == expected_size_fb2);
  REQUIRE(ind_calc2->get_factor_base()->get_size_fb() == expected_size_fb2);
  REQUIRE(ind_calc2->get_factor_base()->get_bound() == expected_bound_fb2);
  REQUIRE(ind_calc2->get_relation_generator()->get_max_num_tests() == expected_max_num_tests2);
}

#endif

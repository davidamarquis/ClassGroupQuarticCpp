#pragma once

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstdint> // for fixed-width int types
#include <cassert>
#include <cmath> // for pow() function
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <antic/nf.h>
#include <antic/nf_elem.h>
#include <flint/arith.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

extern nf_t antic_nf;
//This file is responsible for compiling the source code of Catch and linking together separate test files.
//See here for more information: https://github.com/catchorg/Catch2/blob/v2.x/docs/tutorial.md#scaling-up

// In a Catch project with multiple files, dedicate one file to compile the
// source code of Catch itself and reuse the resulting object file for linking.

// Let Catch provide main():
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

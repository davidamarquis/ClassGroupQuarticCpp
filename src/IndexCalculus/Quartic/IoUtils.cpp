#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <antic/nf.h>
#include <antic/nf_elem.h>
#include <flint/arith.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include "IoUtils.hpp"
#include "include/rapidjson/error/en.h"
#include "include/rapidjson/istreamwrapper.h"
#include "include/rapidjson/stringbuffer.h"
#include "include/rapidjson/writer.h"

using namespace rapidjson;
using namespace std;

void LoadDoc(rapidjson::Document &doc, const std::string filepath) {
  std::ifstream ifs { filepath };
  if ( !ifs.is_open() ) {
    throw std::runtime_error("Could not open file for reading!\n");
  }

  IStreamWrapper isw { ifs };
  doc.ParseStream( isw );

  StringBuffer buffer {};
  Writer<StringBuffer> writer { buffer };
  doc.Accept( writer );

  if ( doc.HasParseError() ) {
    string error = "Error  : " + string(GetParseError_En(doc.GetParseError()))  + "\n" + "Offset : " + to_string(doc.GetErrorOffset()) + '\n';
    throw std::runtime_error(error);
  }
}

void Split(string delimiter, string str, vector<string> &tokens) {
  size_t pos = 0;
  string token;
  while ((pos = str.find(delimiter)) != string::npos) {
    token = str.substr(0, pos);
    tokens.push_back(token);
    str.erase(0, pos + delimiter.length());
  }
  tokens.push_back(str);
}

bool StartsWith(string first, string second) {
  auto res = std::mismatch(first.begin(), first.end(), second.begin());
  return res.first == first.end();
}
#ifndef THESIS_IOUTILS_HPP
#define THESIS_IOUTILS_HPP

#include <string>
#include <vector>
#include "include/rapidjson/document.h"

void LoadDoc(rapidjson::Document &doc, const std::string filepath);
void Split(std::string delimiter, std::string str, std::vector<std::string> &tokens);
bool StartsWith(std::string first, std::string second);

#endif  // THESIS_IOUTILS_HPP
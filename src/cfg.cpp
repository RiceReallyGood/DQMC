#include "cfg.h"

#include <cctype>
#include <fstream>
#include <iostream>

Config::Config(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Can not open file \"" << filename << "\"" << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        // get rid of comments
        line = line.substr(0, line.find('#'));
        // std::cout << line << ";" << std::endl;
        if (line.empty()) continue;
        std::string::size_type idx = line.find(delimiter);
        if (idx == std::string::npos) continue;
        std::string attr = line.substr(0, idx);
        std::string val = line.substr(idx + 1);

        if (ToUpper(strip(attr)) == "" || strip(val) == "") continue;

        // std::cout << attr << ":" << val << ";" << std::endl;
        m[attr] = val;
    }
}

std::string Config::GetCfg(std::string attr) {
    ToUpper(attr);
    if (m.find(attr) == m.end()) return "";
    return m[attr];
}

std::string& Config::ToUpper(std::string& s) {
    for (int i = 0; i < s.length(); i++) {
        s[i] = std::toupper(s[i]);
    }
    return s;
}

std::string& Config::strip(std::string& s) {
    std::string::size_type begIdx, endIdx;
    begIdx = s.find_first_not_of(' ');
    endIdx = s.find_last_not_of(' ');
    s = (begIdx == std::string::npos) ? "" : s.substr(begIdx, endIdx - begIdx + 1);
    return s;
}

Config::~Config() {}
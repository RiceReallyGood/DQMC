#ifndef DQMC_CFG_H
#define DQMC_CFG_H

#include <string>
#include <unordered_map>

class Config {
public:
    Config(const std::string& filename);
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
    ~Config();

    std::string GetCfg(std::string attr);

private:
    std::unordered_map<std::string, std::string> m;
    static std::string& ToUpper(std::string& s);
    static std::string& strip(std::string& s);
    static const char delimiter = '=';
};

#endif
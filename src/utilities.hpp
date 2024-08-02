#pragma once

#include <cstdlib> // for std::exit
#include <iostream>
#include <spdlog/spdlog.h>
#include <string>

inline void exitWithError(const std::string& errorMessage) {
    spdlog::error(errorMessage);
    std::exit(EXIT_FAILURE);
}

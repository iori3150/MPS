#pragma once

#include <cstdlib> // for std::exit
#include <iostream>
#include <string>

inline void exitWithError(const std::string& errorMessage) {
    std::cout << "ERROR: " << errorMessage << std::endl << std::endl;
    std::exit(EXIT_FAILURE);
}

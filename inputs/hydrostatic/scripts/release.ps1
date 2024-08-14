clang++ -std=c++20 -I ../../submodules/eigen -c main.cpp -o build/main.o
clang++ -std=c++20 -I ../../submodules/eigen -c ../../src/particle.cpp -o build/particle.o
clang++ -std=c++20 -I ../../submodules/eigen -I ../../submodules/csv-parser/single_include -I ../../submodules/spdlog/include -c ../../src/exporter.cpp -o build/exporter.o
clang++ build/main.o build/particle.o build/exporter.o -o build/main.exe
clang++ -std=c++20 -I ../eigen-3.4.0 -c main.cpp -o build/main.o
clang++ -std=c++20 -I ../eigen-3.4.0 -c ../src/particle.cpp -o build/particle.o
clang++ -std=c++20 -I ../eigen-3.4.0 -I ../csv-parser-2.3.0/include -c ../src/exporter.cpp -o build/exporter.o
clang++ build/main.o build/particle.o build/exporter.o -o build/main.exe
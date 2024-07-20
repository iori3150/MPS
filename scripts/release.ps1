Remove-Item build/*
clang++ -std=c++17 -fopenmp -I eigen-3.4.0 -c src/main.cpp -o build/main.o
clang++ -std=c++17 -fopenmp -I eigen-3.4.0 -c src/simulation.cpp -o build/simulation.o
clang++ -std=c++17 -fopenmp -I eigen-3.4.0 -c src/particle.cpp -o build/particle.o
clang++ -fopenmp build/main.o build/simulation.o build/particle.o -o build/main.exe
Remove-Item build/*
clang++ -std=c++17 -g -I eigen-3.4.0 -c src/main.cpp -o build/main.o
clang++ -std=c++17 -g -I eigen-3.4.0 -c src/simulation.cpp -o build/simulation.o
clang++ -g build/main.o build/simulation.o -o build/main.exe
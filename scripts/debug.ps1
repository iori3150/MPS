Remove-Item build/*
clang++ -std=c++20 -g -I eigen-3.4.0 -I csv-parser-2.3.0/include -c src/main.cpp -o build/main.o
clang++ -std=c++20 -g -I eigen-3.4.0 -I csv-parser-2.3.0/include -c src/simulation.cpp -o build/simulation.o
clang++ -std=c++20 -g -I eigen-3.4.0 -I csv-parser-2.3.0/include -c src/saver.cpp -o build/saver.o
clang++ -std=c++20 -g -fopenmp -I eigen-3.4.0 -c src/particle.cpp -o build/particle.o
clang++ -std=c++20 -g -fopenmp -I eigen-3.4.0 -c src/mps.cpp -o build/mps.o
clang++ -std=c++20 -g -fopenmp -I eigen-3.4.0 -c src/bucket.cpp -o build/bucket.o
clang++ -fopenmp -g build/main.o build/simulation.o build/particle.o build/mps.o build/bucket.o build/saver.o -o build/main.exe
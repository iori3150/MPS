# Create build folder if it does not exist
if (-Not (Test-Path -Path build)) {
  New-Item -ItemType Directory -Path build
}
Remove-Item build/*
clang++ -std=c++20 -g -I eigen-3.4.0 -I csv-parser-2.3.0/include -c src/main.cpp -o build/main.o
clang++ -std=c++20 -g -I eigen-3.4.0 -I csv-parser-2.3.0/include -c src/simulation.cpp -o build/simulation.o
clang++ -std=c++20 -g -I eigen-3.4.0 -I csv-parser-2.3.0/include -c src/exporter.cpp -o build/exporter.o
clang++ -std=c++20 -g -I eigen-3.4.0 -I csv-parser-2.3.0/include -c src/particle.cpp -o build/particle.o
clang++ -std=c++20 -g -I eigen-3.4.0 -I csv-parser-2.3.0/include -fopenmp -c src/mps.cpp -o build/mps.o
clang++ -std=c++20 -g -I eigen-3.4.0 -I csv-parser-2.3.0/include -c src/bucket.cpp -o build/bucket.o
clang++ -g -fopenmp build/main.o build/simulation.o build/particle.o build/mps.o build/bucket.o build/exporter.o -o build/main.exe
del main\build_error.log
del main\main.exe
clang++ main\main.cpp -std=c++17 -fopenmp -I eigen-3.4.0 -o main\main.exe 2>main\build_error.log
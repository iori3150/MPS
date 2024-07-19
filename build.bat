del main\build_error.log
del main\main.exe
g++ main\main.cpp -std=c++17 -fopenmp -I C:\\eigen-3.4.0 -o main\main.exe 2>main\build_error.log
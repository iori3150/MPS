del result\*.prof
.\main\main.exe 2>result\error.log

del prof2vtu\prof2vtu.exe
g++ prof2vtu\prof2vtu.cpp -std=c++17  -o prof2vtu\prof2vtu.exe

del result\*.vtu
.\prof2vtu\prof2vtu.exe 2>result\error.log
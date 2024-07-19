SET /p answer="Execute main.cpp? (Y/N)"
if /i {%ANSWER%}=={y} (goto :yes)
if /i {%ANSWER%}=={Y} (goto :yes)
if /i {%ANSWER%}=={yes} (goto :yes)
EXIT

:yes
del result\prof\*.prof
del result\vtu\*.vtu
.\main\main.exe 2> result\error.log
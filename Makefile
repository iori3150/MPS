# Define compiler and flags
CXX = clang++
CXXFLAGS = -std=c++20 -O2 -I eigen-3.4.0 -I csv-parser-2.3.0/single_include

# Define directories
SRCDIR = src
BUILDDIR = build

TARGET = $(BUILDDIR)/main.exe

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# All object files
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SOURCES))
# Object files that require OpenMP
OPENMP_OBJECTS = $(BUILDDIR)/mps.o
# Define other object files by filtering out OpenMP object files
OTHER_OBJECTS = $(filter-out $(OPENMP_OBJECTS), $(OBJECTS))

# ---------------------
# Target: all (default)
# ---------------------
all: $(TARGET)

# Link all object files and create main.exe
$(TARGET): $(OBJECTS)
	$(CXX) -fopenmp $(OTHER_OBJECTS) $(OPENMP_OBJECTS) -o $@

# Compile source files that require OpenMP
$(OPENMP_OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@powershell -Command "if (!(Test-Path '$(BUILDDIR)')) { New-Item -ItemType Directory -Path '$(BUILDDIR)' }"
	$(CXX) $(CXXFLAGS) -fopenmp -c $< -o $@

# Compile source files that do not require OpenMP
$(OTHER_OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@powershell -Command "if (!(Test-Path '$(BUILDDIR)')) { New-Item -ItemType Directory -Path '$(BUILDDIR)' }"
	$(CXX) $(CXXFLAGS) -c $< -o $@


# -------------
# Target: clean
# -------------
clean:
	@powershell -Command "if (Test-Path '$(BUILDDIR)') { Remove-Item '$(BUILDDIR)/*' }"

# Declare 'all' and 'clean' as phony targets
.PHONY: all clean

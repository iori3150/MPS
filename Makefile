# How to use ------------------------------------
# make							-> compile in release mode
# make MODE=debug		-> compile in debug mode
# make clean 				-> remove the build directory
# -----------------------------------------------

# compiler, compiler flags, linker flags
CXX = clang++
CXXFLAGS = -std=c++20 -I eigen-3.4.0 -I csv-parser-2.3.0/single_include -I fkYAML-0.3.9/single_include -I spdlog/include
CXXFLAGS_OPENMP = -fopenmp
LDFLAGS = -fopenmp

# build mode: release (default) or debug
MODE ?= release
ifeq ($(MODE),debug)
	CXXFLAGS += -g -O0
else
	CXXFLAGS += -O2
endif

SRCDIR = src
BUILDDIR = build
BUILDMODEDIR = $(BUILDDIR)/$(MODE)
TARGET = $(BUILDMODEDIR)/main.exe
SOURCES = $(wildcard $(SRCDIR)/*.cpp) # Source files to compile
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(BUILDMODEDIR)/%.o,$(SOURCES)) # Object files to link

# ---------------------
# Target: all (default)
# ---------------------
all: $(TARGET)

$(TARGET): $(OBJECTS)
  # Link object files to create executable
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

$(BUILDMODEDIR)/%.o: $(SRCDIR)/%.cpp
  # Create build directory if it doesn't exist
	@powershell -Command "if (!(Test-Path '$(BUILDMODEDIR)')) { New-Item -ItemType Directory -Path '$(BUILDMODEDIR)' }"

ifeq ($(notdir $<),mps.cpp)
  # Compile mps.cpp with OpenMP
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_OPENMP) -c $< -o $@
else
  # Compile other source files into object files
	$(CXX) $(CXXFLAGS) -c $< -o $@
endif

# ---------------------
# Target: clean
# ---------------------
clean:
	@powershell -Command "if (Test-Path '$(BUILDDIR)') { Remove-Item -Recurse -Force '$(BUILDDIR)/*' }"

# Declare 'all' and 'clean' as phony targets
.PHONY: all clean

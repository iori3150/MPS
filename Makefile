CXX = clang++ # compiler
CXXFLAGS = -std=c++20 -O2 -I eigen-3.4.0 -I csv-parser-2.3.0/single_include -fopenmp # compiler flags
LDFLAGS = -fopenmp # linker flags

# Define directories
SRCDIR = src
BUILDDIR = build

TARGET = $(BUILDDIR)/main.exe

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# All object files
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SOURCES))

# ---------------------
# Target: all (default)
# ---------------------
all: $(TARGET)

$(TARGET): $(OBJECTS)

  # Create build directory if it doesn't exist
	@powershell -Command "if (!(Test-Path '$(BUILDDIR)')) { New-Item -ItemType Directory -Path '$(BUILDDIR)' }"

  # Link object files to create executable
	$(CXX) -fopenmp $(OBJECTS) -o $@

# Compile each source file into an object file
$(OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# -------------
# Target: clean
# -------------
clean:
	@powershell -Command "if (Test-Path '$(BUILDDIR)') { Remove-Item '$(BUILDDIR)/*' }"

# Declare 'all' and 'clean' as phony targets
.PHONY: all clean

# How to use ------------------------------------
# make							-> compile in release mode
# make MODE=debug		-> compile in debug mode
# make clean 				-> remove the build directory
# -----------------------------------------------

# compiler, compiler flags, linker flags
CXX = clang++
CXXFLAGS = -std=c++20 -I eigen-3.4.0 -I csv-parser-2.3.0/single_include -fopenmp
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
  # Create build directory if it doesn't exist
	@powershell -Command "if (!(Test-Path '$(BUILDMODEDIR)')) { New-Item -ItemType Directory -Path '$(BUILDMODEDIR)' }"

  # Link object files to create executable
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

# Compile each source file into an object file
$(BUILDMODEDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---------------------
# Target: clean
# ---------------------
clean:
	@powershell -Command "if (Test-Path '$(BUILDDIR)') { Remove-Item -Recurse -Force '$(BUILDDIR)/*' }"

# Declare 'all' and 'clean' as phony targets
.PHONY: all clean

# How to use ------------------------------------
# make							-> compile in release mode
# make MODE=debug		-> compile in debug mode
# make clean 				-> remove the build directory
# -----------------------------------------------

CXX = clang++
CXXFLAGS = -std=c++20 -I submodules/eigen -I submodules/csv-parser/single_include -I submodules/fkYAML/single_include -I submodules/spdlog/include
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
  # Mac or Linux
	@if [ "$(shell uname)" = "Darwin" ] || [ "$(shell uname)" = "Linux" ]; then \
		mkdir -p $(BUILDMODEDIR); \

  # Windows Command Prompt
	elif [ -n "$(ComSpec)" ]; then \ 
	  if not exist "$(BUILDMODEDIR)" mkdir "$(BUILDMODEDIR)"; \

  # Windows PowerShell
	else \
		powershell -Command "if (!(Test-Path '$(BUILDMODEDIR)')) { New-Item -ItemType Directory -Path '$(BUILDMODEDIR)' }"; \
	fi

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
  # Mac or Linux
	@if [ "$(shell uname)" = "Darwin" ] || [ "$(shell uname)" = "Linux" ]; then \
		rm -rf $(BUILDDIR); \

  # Windows Command Prompt
	elif [ -n "$(ComSpec)" ]; then \
		if exist "$(BUILDDIR)" rmdir /s /q "$(BUILDDIR)"; \

  # Windows PowerShell
  else \
		powershell -Command "if (Test-Path '$(BUILDDIR)') { Remove-Item -Recurse -Force '$(BUILDDIR)/*' }"; \
  fi

# Declare 'all' and 'clean' as phony targets
.PHONY: all clean

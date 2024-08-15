CXX = clang++
INCLUDES = \
	-I ../../submodules/eigen \
	-I ../../submodules/spdlog/include \
	-I ../../submodules/csv-parser/single_include
CXXFLAGS = -std=c++20 $(INCLUDES)

BUILDDIR = build
TARGET = $(BUILDDIR)/main.exe
SOURCES = main.cpp ../../src/particle.cpp ../../src/exporter.cpp
OBJECTS = $(patsubst %.cpp,$(BUILDDIR)/%.o,$(notdir $(SOURCES)))

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@

$(BUILDDIR)/%.o: ../../src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR)/main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: all

CXX = g++
CXXFLAGS = -g -O2 -I./pcm/src/ -DPCM_DYNAMIC_LIB -std=c++17 -fopenmp # Add -std=c++17 or your desired standard
LDFLAGS = -L./pcm/build/lib/ -Wl,-rpath=./pcm/build/lib/ -fopenmp
LIBS = -lpcm -ldl

SOURCE ?= testes.cpp
TARGET = $(basename $(notdir $(SOURCE)))  # Use filename without extension
SOURCES = $(SOURCE) 
OBJECTS = $(SOURCES:.cpp=.o)

# Default rule:  What to build when you just type "make"
$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $(TARGET)
	rm -f $(OBJECTS)
	
# Rule to compile .cpp files to .o files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

#  Cleanup rule
clean:
	rm -f $(TARGET) $(OBJECTS)

.PHONY: clean #  Declare clean as not being a real file.

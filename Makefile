CXX=/opt/homebrew/opt/llvm/bin/clang++

CPPFLAGS = --std=c++14 -fopenmp

CXXFLAGS = -std=c++14 -O2 -fopenmp \
           -I./include \
           -I/Users/hous/Github/eigen \
           -I/opt/homebrew/opt/libomp/include

LDFLAGS = -L/opt/homebrew/opt/libomp/lib
LDLIBS = -lomp

SRC = $(wildcard *.cpp)
OUT = main.out


all: $(OUT)
$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT) $(LDFLAGS) $(LDLIBS)

# profile: $(OUT)
# $(OUT): $(SRC)
# 	xcrun $(CXX) $(CXXFLAGS) $(SRC) -o $(OUT) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OUT)



# SOURCES := $(wildcard *.cpp)
# OBJECTS := $(SOURCES:.cpp=.o)

# BASE := simple

# HEADERS := $(wildcard *.H)

# PYBIND_INCLUDES := $(shell python -m pybind11 --includes)

# PYBIND_SUFFIX := $(shell python3-config --extension-suffix)

# PYTHON_LIBRARY := ${BASE}${PYBIND_SUFFIX}

# ALL: ${PYTHON_LIBRARY}

# CXXFLAGS := -O3  -Wall -Wextra -shared -std=c++17 -fPIC ${PYBIND_INCLUDES}

# %.o : %.cpp
# 	g++ ${CXXFLAGS} -c $<

# ${PYTHON_LIBRARY}: ${OBJECTS} ${HEADERS}
# 	$(CXX) -O3  -Wall -Wextra -shared -std=c++17 -fPIC ${PYBIND_INCLUDES} simple.cpp -o $@ $<

# print-%: ; @echo $* is $($*)
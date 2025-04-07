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
	xcrun $(CXX) $(CXXFLAGS) $(SRC) -o $(OUT) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OUT)

CXX=/opt/homebrew/opt/llvm/bin/clang++

CXXFLAGS = -std=c++14 -O2 -fopenmp \
           -I./include \
           -I../extern \
           -I/opt/homebrew/opt/libomp/include \
		   -I

LDFLAGS = -L/opt/homebrew/opt/libomp/lib
LDLIBS = -lomp

SRC = $(wildcard *.cpp)
OUT = main.out

all:
	$(CXX) $(CXXFLAGS) -g emd_wrap.cpp solver.cpp solver_memory.cpp utils.cpp main.cpp -o main.out $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OUT)
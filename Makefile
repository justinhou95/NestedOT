main.out: main.cpp solver.cpp utils.cpp generate_path.cpp
	g++ --std=c++14 main.cpp solver.cpp utils.cpp generate_path.cpp -o main.out
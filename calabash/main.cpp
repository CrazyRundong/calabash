#include "calabash.h"
#include <iostream>
#include <fstream>
#include <random>

int main(int argc, char** argv) {
#ifdef _DEBUG_DIAMOND_
	// test edges
	std::ifstream ifs3("../input/2.txt");
	// Edges eg3(ifs3);

	// test calcbash
	// Diamond dia3(eg3, eg3.num_nodes());
	Diamond dia3(ifs3);
	std::cout << "Diamond score:\t" << dia3.get_power() << std::endl
		<< "now randomly change 50 states, new score:\t";
	std::uniform_int_distribution<int> u(1, dia3.get_num_nodes());
	std::default_random_engine re;
	for (int i = 0; i < 50; ++i) {
		dia3.set_state(u(re));
	}
	std::cout << dia3.get_power() << std::endl;

	// test memory when copy
	Diamond dia3cp(dia3);
	for (int i = 0; i < 50; ++i) {
		dia3cp.set_state(u(re));
	}
	std::cout << "After copy:" << std::endl
		<< "Dia3   score:\t" << dia3.get_power() << std::endl
		<< "Dia3cp score:\t" << dia3cp.get_power() << std::endl;
#endif

	int numStep = std::stoi(argv[1]);
	std::string graphFilePath(argv[2]);
	std::string resultFilePath(argv[3]);

	RandomWalkingSolver solver(graphFilePath);
	std::ofstream resultFile(resultFilePath, std::ofstream::out);

	std::cout << "Diamond power:\t" << solver.solve(numStep) << std::endl;
	solver.print_result(resultFile);

	return 0;
}

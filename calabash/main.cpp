#include "calabash.h"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

int main(int argc, char** argv) {
	// test edges
	std::ifstream ifs3("../input/3.txt");
	Edges eg3(ifs3);

	// test calcbash
	Diamond dia3(eg3, eg3.num_nodes());
	std::cout << "Diamond score:\t" << dia3.get_power() << std::endl;

	return 0;
}
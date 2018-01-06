#include "calabash.h"
#include <iostream>
#include <fstream>
#include <random>

int main(int argc, char** argv) {
    double bestPower;
    int numStep = std::stoi(argv[1]);
    std::string graphFilePath(argv[2]);
    std::string resultFilePath(argv[3]);

    RandomWalkingSolver solver(graphFilePath);
    bestPower = solver.solve(numStep);

#ifdef _DEBUG_DIAMOND_
    std::cout << "Diamond power:\t" << bestPower << std::endl;
    std::ofstream resultFile(resultFilePath, std::ofstream::out);
    solver.print_result(resultFile);
#else
    solver.print_result(std::cout);
#endif // _DEBUG_DIAMOND_

    return 0;
}

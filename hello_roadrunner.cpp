#include <iostream>
#include "rr/rrRoadRunner.h"

int main() {
    std::cout << "Hello, World!" << std::endl;
    rr::RoadRunner rr("./glycolysis.xml");
    return 0;
}
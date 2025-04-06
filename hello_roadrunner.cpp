#include <iostream>
#include <filesystem>
#include <fstream>
#include "rr/rrRoadRunner.h"
#include <sbml/SBMLReader.h>
#include <sbml/SBMLDocument.h>

int main() {
    const std::string path = "glycolysis.xml";

    std::cout << "Checking for file: " << path << std::endl;
    std::cout << "Working directory: " << std::filesystem::current_path() << std::endl;

    if (!std::filesystem::exists(path)) {
        std::cerr << "ERROR: File not found!" << std::endl;
        return 1;
    }

    std::ifstream file("glycolysis.xml");
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string sbml = buffer.str();

    auto document = libsbml::readSBMLFromString(sbml.c_str());
    if (document->getNumErrors() > 0) {
        document->printErrors();  // Shows all parse/validation errors
        return 1;
    }


    try {
        rr::RoadRunner rr;
        rr.load(path);
        std::cout << "Model loaded successfully!" << std::endl;
        rr::SimulateOptions& opts = rr.getSimulateOptions();
        opts.start = 0;
        opts.duration = 100;
        opts.steps = 100;
        const ls::DoubleMatrix* result = rr.simulate();

        if (result) {
          std::cout << "Simulation completed. First 5 rows:" << std::endl;
          for (int i = 0; i < std::min<size_t>(5, result->numRows()); ++i) {
              for (int j = 0; j < result->numCols(); ++j) {
                  std::cout << (*result)(i, j) << '\t';
              }
              std::cout << '\n';
          }
      } else {
          std::cerr << "Simulation failed or returned no data." << std::endl;
      }
      

    } catch (const std::exception& e) {
        std::cerr << "Error loading model: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

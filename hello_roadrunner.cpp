#include <iostream>
#include <conio.h>
#include <iomanip> // required for setw, setprecision, etc.
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


    auto str = rr::RoadRunner::getExtendedVersionInfo();
    std::cout << "Version ptr: " << static_cast<const void*>(str.c_str()) << std::endl;
    std::cout << "LibRoadRunner version: " << rr::RoadRunner::getExtendedVersionInfo() << std::endl;


    try {
        rr::RoadRunner rr;
        // rr.setIntegrator("rk4");
        rr.load(path);
        std::cout << "Model loaded successfully!" << std::endl;

        std::cout << "Integrator: " << rr.getIntegrator() << std::endl;
        std::cout << "Number of reactions: " << rr.getNumberOfReactions() << std::endl;
        std::cout << "Number of compartments: " << rr.getNumberOfCompartments() << std::endl;
        std::cout << "Number of floating species: " << rr.getNumberOfFloatingSpecies() << std::endl;
        std::cout << "Number of boundary species: " << rr.getNumberOfBoundarySpecies() << std::endl;

        auto ids = rr.getFloatingSpeciesIds();
        auto concs = rr.getFloatingSpeciesConcentrationsV();
        for (size_t i = 0; i < ids.size(); ++i) {
            std::cout << ids[i] << ": " << concs[i] << std::endl;
        }

        for (const auto& id : ids) {
            std::cout << id << ": " << rr.getValue(id) << std::endl;
        }

        auto boundaryIds = rr.getBoundarySpeciesIds();
        for (const auto& id : boundaryIds) {
            std::cout << "(boundary) " << id << ": " << rr.getValue(id) << std::endl;
        }

        if (!rr.getValue("P") || !rr.getValue("ATP") || !rr.getValue("ADP")) {
            std::cerr << "ERROR: One or more selected species do not exist." << std::endl;
            return 1;
        }

        std::cout << "Running real-time simulation. Press 'q' to quit.\n";

        double t = 0.0;
        double dt = 0.001;

        // Set useful output variables
        rr.setSelections({"time", "GLCi", "ATP", "ADP", "P", "NAD", "NADH", "PYR"}); // Replace with actual species IDs if needed

        // std::stringstream* state = rr.saveStateS();  // save initial state

        while (true) {
            // Check for keypress
            if (_kbhit()) {
                char ch = _getch();
                std::cout << "\nPressed: " << ch << std::endl;
                if (ch == 'q' || ch == 'Q') {
                    std::cout << "Exiting simulation.\n";
                    break;
                } else if (ch == 'p' || ch == 'P') {
                    rr.oneStep(t, dt); // advance
                    t += dt;
                } else if (ch == '1') {
                    // rr.loadStateS(state);
                    double glci = rr.getValue("GLCi");
                    // rr.reset(rr::SelectionRecord::ALL);
                    // rr.getIntegrator()->restart(t);
                    rr.setValue("GLCi", glci + 1.0);
                    std::cout << "GLCi increased to " << (glci + 1.0) << std::endl;
                } else if (ch == 'i' || ch == 'I') {
                    auto reactionRates = rr.getReactionRates();
                    auto reactionIds = rr.getReactionIds();
                    std::cout << "Reaction rates:\n";
                    for (size_t i = 0; i < reactionRates.size(); ++i) {
                        std::cout << reactionIds[i] << ": " << reactionRates[i] << std::endl;
                    }
                } else if (ch == 'r' || ch == 'R') {
                    // rr.loadStateS(state);
                    rr.reset(rr::SelectionRecord::ALL); // reset all
                }
                auto values = rr.getSelectedValues();
    
                std::cout << std::fixed << std::setprecision(2); // fixed-point notation with 2 decimal places
                std::cout << "t=" << std::setw(5) << t << " [ ";
                for (double v : values) {
                    std::cout << std::setw(8) << v << " ";
                }
                std::cout << "]" << std::endl;
    
                // Optional: sleep to simulate real time
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            } else {    
                // rr.loadStateS(state);
            }
            // state = nullptr;


            // Save new state for next frame
            // state = rr.saveStateS();

            // rr.oneStep(t, dt);
            // t += dt;

        }

        // rr::SimulateOptions& opts = rr.getSimulateOptions();
        // opts.start = 0;
        // opts.duration = 100;
        // opts.steps = 100;
        // const ls::DoubleMatrix* result = rr.simulate();

        // if (result) {
        //   std::cout << "Simulation completed. First 5 rows:" << std::endl;
        //   for (int i = 0; i < std::min<size_t>(5, result->numRows()); ++i) {
        //       for (int j = 0; j < result->numCols(); ++j) {
        //           std::cout << (*result)(i, j) << '\t';
        //       }
        //       std::cout << '\n';
        //   }
        // } else {
        //     std::cerr << "Simulation failed or returned no data." << std::endl;
        // }
      

    } catch (const std::exception& e) {
        std::cerr << "Error loading model: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

#include <iostream>
#include <conio.h>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include "rr/rrRoadRunner.h"
#include <sbml/SBMLReader.h>
#include <sbml/SBMLDocument.h>
#include <map> // For tracking previous values

// Function to check if the file exists
bool checkFileExists(const std::string& path) {
    std::cout << "Checking for file: " << path << std::endl;
    std::cout << "Working directory: " << std::filesystem::current_path() << std::endl;

    if (!std::filesystem::exists(path)) {
        std::cerr << "ERROR: File not found!" << std::endl;
        return false;
    }
    return true;
}

// Function to load SBML content from a file
std::string loadSBMLFromFile(const std::string& path) {
    std::ifstream file(path);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

// Function to validate SBML content
bool validateSBML(const std::string& sbml) {
    auto document = libsbml::readSBMLFromString(sbml.c_str());
    if (document->getNumErrors() > 0) {
        document->printErrors();
        return false;
    }
    return true;
}

// Function to log all species and their concentrations with color coding
void logSpeciesWithColor(rr::RoadRunner& rr, std::map<std::string, double>& previousValues) {
    auto ids = rr.getFloatingSpeciesIds();
    auto concs = rr.getFloatingSpeciesConcentrationsV();

    std::cout << "Floating species concentrations:\n";
    for (size_t i = 0; i < ids.size(); ++i) {
        const std::string& id = ids[i];
        double currentValue = concs[i];
        double previousValue = previousValues[id];
        double diff = currentValue - previousValue;

        if (diff > 0) {
            std::cout << "\033[32m"; // Green for increase
        } else if (diff < 0) {
            std::cout << "\033[31m"; // Red for decrease
        } else {
            std::cout << "\033[0m"; // Reset color
        }

        std::cout << id << ": " << currentValue << "(" << diff << ")" << "\033[0m" << std::endl; // Reset color after value
        previousValues[id] = currentValue; // Update previous value
    }

    auto boundaryIds = rr.getBoundarySpeciesIds();
    std::cout << "Boundary species concentrations:\n";
    for (const auto& id : boundaryIds) {
        double currentValue = rr.getValue(id);
        double previousValue = previousValues[id];

        if (currentValue > previousValue) {
            std::cout << "\033[32m"; // Green for increase
        } else if (currentValue < previousValue) {
            std::cout << "\033[31m"; // Red for decrease
        } else {
            std::cout << "\033[0m"; // Reset color
        }

        std::cout << id << ": " << currentValue << "\033[0m" << std::endl; // Reset color after value
        previousValues[id] = currentValue; // Update previous value
    }
}

// Function to increase the concentration of a chosen species
void increaseSpecies(rr::RoadRunner& rr, const std::string& speciesId, double increment) {
    try {
        double currentValue = rr.getValue(speciesId);
        rr.setValue(speciesId, currentValue + increment);
        std::cout << speciesId << " increased to " << (currentValue + increment) << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to modify species " << speciesId << ": " << e.what() << std::endl;
    }
}

// Function to merge two SBML models
std::string mergeSBMLModels(const std::string& sbml1, const std::string& sbml2) {
    auto doc1 = libsbml::readSBMLFromString(sbml1.c_str());
    auto doc2 = libsbml::readSBMLFromString(sbml2.c_str());

    if (!doc1 || !doc2) {
        throw std::runtime_error("Error reading SBML documents.");
    }

    auto model1 = doc1->getModel();
    auto model2 = doc2->getModel();

    if (!model1 || !model2) {
        throw std::runtime_error("Error retrieving models from SBML documents.");
    }

    // Rename compartment "cytoplasm" to "cytosol" in both models
    auto renameCompartment = [](libsbml::Model* model) {
        auto compartment = model->getCompartment("cytoplasm");
        if (compartment) {
            compartment->setId("cytosol");
        }
        for (unsigned int i = 0; i < model->getNumSpecies(); ++i) {
            auto species = model->getSpecies(i);
            if (species->getCompartment() == "cytoplasm") {
                species->setCompartment("cytosol");
            }
        }
    };
    renameCompartment(model1);
    renameCompartment(model2);

    // Merge compartments from model2 into model1
    for (unsigned int i = 0; i < model2->getNumCompartments(); ++i) {
        auto compartment = model2->getCompartment(i)->clone();
        if (!model1->getCompartment(compartment->getIdAttribute())) {
            model1->addCompartment(static_cast<libsbml::Compartment*>(compartment));
        } else {
            delete compartment; // Avoid memory leak
        }
    }

    // Remove "cytoplasm" if it is unused in model1
    if (model1->getCompartment("cytoplasm")) {
        bool isUsed = false;
        for (unsigned int i = 0; i < model1->getNumSpecies(); ++i) {
            if (model1->getSpecies(i)->getCompartment() == "cytoplasm") {
                isUsed = true;
                break;
            }
        }
        if (!isUsed) {
            model1->removeCompartment("cytoplasm");
        }
    }

    // Merge species from model2 into model1
    for (unsigned int i = 0; i < model2->getNumSpecies(); ++i) {
        auto species = model2->getSpecies(i)->clone();
        if (!model1->getSpecies(species->getIdAttribute())) {
            model1->addSpecies(static_cast<libsbml::Species*>(species));
        } else {
            delete species; // Avoid memory leak
        }
    }

    // Merge reactions from model2 into model1
    for (unsigned int i = 0; i < model2->getNumReactions(); ++i) {
        auto reaction = model2->getReaction(i)->clone();
        if (!model1->getReaction(reaction->getIdAttribute())) {
            model1->addReaction(static_cast<libsbml::Reaction*>(reaction));
        } else {
            delete reaction; // Avoid memory leak
        }
    }

    // Merge parameters from model2 into model1
    for (unsigned int i = 0; i < model2->getNumParameters(); ++i) {
        auto parameter = model2->getParameter(i)->clone();
        if (!model1->getParameter(parameter->getIdAttribute())) {
            model1->addParameter(static_cast<libsbml::Parameter*>(parameter));
        } else {
            delete parameter; // Avoid memory leak
        }
    }

    // Ensure all parameters in model1 have initial values
    for (unsigned int i = 0; i < model1->getNumParameters(); ++i) {
        auto parameter = model1->getParameter(i);
        if (!parameter->isSetValue()) { // Check only if the value is set
            parameter->setValue(0.0); // Assign a default value of 0.0
        }
    }

    // Write the merged model to a string
    libsbml::SBMLWriter writer;
    std::stringstream mergedSBML;
    writer.writeSBML(doc1, mergedSBML);

    delete doc1;
    delete doc2;

    return mergedSBML.str();
}

int main() {
    const std::string glyPath = "glycolysis.xml";
    const std::string tcaPath = "tca_cycle.xml";

    if (!checkFileExists(glyPath) || !checkFileExists(tcaPath)) {
        return 1;
    }

    try {
        // Load and merge SBML models
        std::string sbml1 = loadSBMLFromFile(glyPath);
        std::string sbml2 = loadSBMLFromFile(tcaPath);
        std::string mergedSBML = mergeSBMLModels(sbml1, sbml2);

        // Save the merged SBML to a temporary file
        const std::string mergedPath = "merged_model.xml";
        std::ofstream mergedFile(mergedPath);
        mergedFile << mergedSBML;
        mergedFile.close();

        if (!validateSBML(mergedSBML)) {
            return 1;
        }

        rr::RoadRunner rr;
        rr.load(mergedPath);
        std::cout << "Merged model loaded successfully!" << std::endl;

        std::cout << "Running real-time simulation. Press 'q' to quit.\n";

        double t = 0.0;
        double dt = 0.01;

        // Dynamically set selections to include all floating species
        auto floatingSpeciesIds = rr.getFloatingSpeciesIds();
        std::vector<std::string> selections = {"time"};
        selections.insert(selections.end(), floatingSpeciesIds.begin(), floatingSpeciesIds.end());
        rr.setSelections(selections);

        std::map<std::string, double> previousValues; // Map to track previous values
        for (const auto& id : floatingSpeciesIds) {
            previousValues[id] = rr.getValue(id); // Initialize with current values
        }
        for (const auto& id : rr.getBoundarySpeciesIds()) {
            previousValues[id] = rr.getValue(id); // Initialize with current values
        }

        while (true) {
            if (_kbhit()) {
                char ch = _getch();
                std::cout << "\nPressed: " << ch << std::endl;
                if (ch == 'q' || ch == 'Q') {
                    std::cout << "Exiting simulation.\n";
                    break;
                } else if (ch == 'p' || ch == 'P') {
                    rr.oneStep(t, dt);
                    logSpeciesWithColor(rr, previousValues);
                    t += dt;
                } else if (ch == 'r' || ch == 'R') {
                    rr.reset(rr::SelectionRecord::ALL);
                } else if (ch == 's' || ch == 'S') {
                    std::cout << "Enter species ID to increase: ";
                    std::string speciesId;
                    std::cin >> speciesId;
                    std::cout << "Enter increment value: ";
                    double increment;
                    std::cin >> increment;
                    increaseSpecies(rr, speciesId, increment);
                }

                // auto values = rr.getSelectedValues();
                // std::cout << std::fixed << std::setprecision(2);
                // std::cout << "t=" << std::setw(5) << t << " [ ";
                // for (double v : values) {
                //     std::cout << std::setw(8) << v << " ";
                // }
                // std::cout << "]" << std::endl;

                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

#include <iostream>
#include <conio.h>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include "rr/rrRoadRunner.h"
#include <sbml/SBMLReader.h>
#include <sbml/SBMLDocument.h>
#include <map> // For tracking previous values
#include <numeric> // For std::iota
#include <sbml/SBMLWriter.h>
#include <sbml/SBMLReader.h>
#include <sbml/SBMLDocument.h>
#include <sbml/Model.h>
#include <rr/rrExecutableModel.h>

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

    // Access the ExecutableModel from RoadRunner
    auto* executableModel = rr.getModel();
    if (!executableModel) {
        throw std::runtime_error("No model loaded in RoadRunner.");
    }

    std::cout << "Floating species concentrations:\n";
    for (size_t i = 0; i < ids.size(); ++i) {
        const std::string& id = ids[i];
        double currentValue = concs[i];
        double previousValue = previousValues[id];
        double diff = currentValue - previousValue;

        // Retrieve the compartment name for the species
        int speciesIndex = executableModel->getFloatingSpeciesIndex(id);
        int compartmentIndex = executableModel->getCompartmentIndexForFloatingSpecies(speciesIndex);
        std::string compartment = executableModel->getCompartmentId(compartmentIndex);

        if (std::abs(diff) >= 0.00005) {
            if (diff > 0) {
                std::cout << "\033[32m"; // Green for increase
            } else {
                std::cout << "\033[31m"; // Red for decrease
            }
        } else {
            std::cout << "\033[0m"; // Reset color
        }

        std::cout << id << " (" << compartment << "): " << std::fixed << std::setprecision(4) << currentValue
                  << " (" << std::fixed << std::setprecision(4) << diff << ")\033[0m" << std::endl; // Reset color after value
        previousValues[id] = currentValue; // Update previous value
    }

    auto boundaryIds = rr.getBoundarySpeciesIds();
    std::cout << "Boundary species concentrations:\n";
    for (const auto& id : boundaryIds) {
        double currentValue = rr.getValue(id);
        double previousValue = previousValues[id];
        double diff = currentValue - previousValue;

        // Retrieve the compartment name for the boundary species
        int speciesIndex = executableModel->getBoundarySpeciesIndex(id);
        int compartmentIndex = executableModel->getCompartmentIndexForBoundarySpecies(speciesIndex);
        std::string compartment = executableModel->getCompartmentId(compartmentIndex);

        if (std::abs(diff) >= 0.00005) {
            if (diff > 0) {
                std::cout << "\033[32m"; // Green for increase
            } else {
                std::cout << "\033[31m"; // Red for decrease
            }
        } else {
            std::cout << "\033[0m"; // Reset color
        }

        std::cout << id << " (" << compartment << "): " << std::fixed << std::setprecision(4) << currentValue << "\033[0m" << std::endl; // Reset color after value
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

// Function to clean up SBML content using libSBML
std::string cleanSBML(const std::string& sbml) {
    auto document = libsbml::readSBMLFromString(sbml.c_str());
    if (!document || document->getNumErrors() > 0) {
        throw std::runtime_error("Error reading SBML document.");
    }

    auto model = document->getModel();
    if (!model) {
        throw std::runtime_error("No model found in SBML document.");
    }

    // Remove <annotation> and <notes> from the model and its components
    model->unsetAnnotation();
    model->unsetNotes();

    for (unsigned int i = 0; i < model->getNumCompartments(); ++i) {
        model->getCompartment(i)->unsetAnnotation();
        model->getCompartment(i)->unsetNotes();
    }

    for (unsigned int i = 0; i < model->getNumSpecies(); ++i) {
        auto species = model->getSpecies(i);
        species->unsetAnnotation();
        species->unsetNotes();
        species->unsetMetaId();
        species->unsetName();
        species->unsetSBOTerm();
    }

    for (unsigned int i = 0; i < model->getNumParameters(); ++i) {
        auto parameter = model->getParameter(i);
        parameter->unsetAnnotation();
        parameter->unsetNotes();
        parameter->unsetMetaId();
        parameter->unsetName();
        parameter->unsetSBOTerm();
    }

    for (unsigned int i = 0; i < model->getNumReactions(); ++i) {
        auto reaction = model->getReaction(i);
        reaction->unsetAnnotation();
        reaction->unsetNotes();
        reaction->unsetMetaId();
        reaction->unsetName();
        reaction->unsetSBOTerm();
    }

    for (unsigned int i = 0; i < model->getNumRules(); ++i) {
        model->getRule(i)->unsetAnnotation();
        model->getRule(i)->unsetNotes();
    }

    for (unsigned int i = 0; i < model->getNumInitialAssignments(); ++i) {
        model->getInitialAssignment(i)->unsetAnnotation();
        model->getInitialAssignment(i)->unsetNotes();
    }

    for (unsigned int i = 0; i < model->getNumEvents(); ++i) {
        auto event = model->getEvent(i);
        event->unsetAnnotation();
        event->unsetNotes();
        event->unsetMetaId();
        event->unsetName();
        event->unsetSBOTerm();
    }

    // Write the cleaned SBML to a string
    libsbml::SBMLWriter writer;
    std::stringstream cleanedSBML;
    writer.writeSBML(document, cleanedSBML);

    delete document;
    return cleanedSBML.str();
}

void clampNegativeConcentrations(rr::RoadRunner& rr) {
    if (auto model = rr.getModel()) {
        int numSpecies = model->getNumFloatingSpecies();

        // Get all assignment rule target IDs
        std::list<std::string> assignmentRuleIds;
        model->getAssignmentRuleIds(assignmentRuleIds);

        // Prepare indices of species NOT governed by assignment rules
        std::vector<int> independentIndices;

        for (int i = 0; i < numSpecies; ++i) {
            std::string speciesId = model->getFloatingSpeciesId(i);

            if (std::find(assignmentRuleIds.begin(), assignmentRuleIds.end(), speciesId) == assignmentRuleIds.end()) {
                // Species is not governed by an assignment rule
                independentIndices.push_back(i);
            }
        }

        // Get and clamp concentrations only for independent species
        std::vector<double> concentrations(independentIndices.size());
        model->getFloatingSpeciesConcentrations(independentIndices.size(), independentIndices.data(), concentrations.data());

        for (auto& conc : concentrations) {
            if (conc < 0.0) {
                conc = 0.00000001;
            }
        }

        model->setFloatingSpeciesConcentrations(independentIndices.size(), independentIndices.data(), concentrations.data());
    }
}

// Function to toggle reaction rates
void toggleReactionRate(rr::RoadRunner& rr, const std::string& reactionId) {
    try {
        // Construct the toggle parameter ID
        std::string toggleParamId = "toggle_" + reactionId;

        // Get the current value of the toggle parameter
        double currentValue = rr.getValue(toggleParamId);

        // Flip the parameter value between 0 and 1
        double newValue = (currentValue == 0.0) ? 1.0 : 0.0;
        rr.setValue(toggleParamId, newValue);

        std::cout << "Parameter " << toggleParamId << " set to " << newValue
                  << " (" << (newValue == 1.0 ? "reaction enabled" : "reaction disabled") << ")." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to toggle parameter for reaction " << reactionId << ": " << e.what() << std::endl;
    }
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

        // Clean the merged SBML
        std::string cleanedSBML = cleanSBML(mergedSBML);

        // Save the cleaned SBML to a temporary file
        const std::string cleanedPath = "cleaned_model.xml";
        std::ofstream cleanedFile(cleanedPath);
        cleanedFile << cleanedSBML;
        cleanedFile.close();

        if (!validateSBML(cleanedSBML)) {
            return 1;
        }

        rr::RoadRunner rr;
        rr.load(cleanedPath);
        std::cout << "Cleaned model loaded successfully!" << std::endl;

        // auto integrator = rr.getIntegrator();
        // Reference below for integrator settings:
        // relative_tolerance
        // absolute_tolerance
        // stiff
        // maximum_bdf_order
        // maximum_adams_order
        // maximum_num_steps
        // maximum_time_step
        // minimum_time_step
        // initial_time_step
        // multiple_steps
        // variable_step_size
        // max_output_rows

        // integrator->setValue("absolute_tolerance", 1e-8);
        // integrator->setValue("relative_tolerance", 1e-6);
        // integrator->setValue("maximum_bdf_order", 5); // Optional but safe
        // integrator->setValue("maximum_adams_order", 12); // Optional

        // integrator->setValue("stiff", true);  // uses BDF
        // integrator->setValue("maximum_num_steps", 10000);  // default is 500
        // integrator->setValue("minimum_step_size", 1e-12);



        std::cout << "Running real-time simulation. Press 'p' to play/pause, 'q' to quit.\n";

        double t = 0.0;
        double dt = 0.01;
        bool isPlaying = false; // Play/pause toggle

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
                    isPlaying = !isPlaying; // Toggle play/pause
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
                } else if (ch == 't' || ch == 'T') {
                    // std::cout << "Enter reaction ID to toggle: ";
                    std::string reactionId = "vAK";
                    // std::cin >> reactionId;
                    toggleReactionRate(rr, reactionId);
                }
            }

            if (isPlaying) {
                try {
                    rr.oneStep(t, dt);
                    t += dt;
                    std::cout << "\033[2J\033[H"; // ANSI escape codes to clear screen and move cursor to top
                } catch (const std::exception& e) {
                    std::cerr << "⚠️ CVODE failed at t = " << t
                              << " with dt = " << dt << "\n"
                              << "Error: " << e.what() << "\n";
                
                    // Try a smaller step size next time, or just continue
                    dt *= 0.5;  // or clamp to some min value
                    if (dt < 1e-10) dt = 0.01;  // reset to default
                    isPlaying = !isPlaying; // Toggle play/pause
                }

                // Clamp all floating species concentrations to zero if they go negative
                clampNegativeConcentrations(rr);

                // Clear the console and log species concentrations
                logSpeciesWithColor(rr, previousValues);

                std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Adjust for desired update rate
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

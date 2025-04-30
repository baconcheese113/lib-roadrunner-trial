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
#include <cmath> // For std::abs
#include <chrono> // For performance timing

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

std::pair<double, double> computeOsmoticPressure(rr::RoadRunner& rr) {
    auto speciesIds = rr.getFloatingSpeciesIds();
    auto concentrations = rr.getFloatingSpeciesConcentrationsV();

    double osmoticPressureCytosol = 0.0;
    double osmoticPressureMito = 0.0;
    std::vector<std::string> includedCytosolSpecies;
    std::vector<std::string> includedMitoSpecies;

    std::set<std::string> excludedSpecies = {"CO2", "He", "ETOH", "SUM_P"};

    auto* executableModel = rr.getModel();
    if (!executableModel) {
        throw std::runtime_error("No model loaded in RoadRunner.");
    }

    for (size_t i = 0; i < speciesIds.size(); ++i) {
        const std::string& id = speciesIds[i];
        if (excludedSpecies.find(id) == excludedSpecies.end()) {
            int speciesIndex = executableModel->getFloatingSpeciesIndex(id);
            int compartmentIndex = executableModel->getCompartmentIndexForFloatingSpecies(speciesIndex);
            std::string compartment = executableModel->getCompartmentId(compartmentIndex);

            if (compartment == "cytosol") {
                osmoticPressureCytosol += concentrations[i];
                includedCytosolSpecies.push_back(id);
            }
            else if (compartment == "mitochondrion") {
                osmoticPressureMito += concentrations[i];
                includedMitoSpecies.push_back(id);
            }
        }
    }

    // Logging species lists
    if (!includedCytosolSpecies.empty()) {
        std::cout << "\033[34mCytosol species contributing to osmotic pressure: ";
        for (size_t i = 0; i < includedCytosolSpecies.size(); ++i) {
            std::cout << includedCytosolSpecies[i];
            if (i < includedCytosolSpecies.size() - 1) std::cout << ", ";
        }
        std::cout << "\033[0m\n";
    }
    if (!includedMitoSpecies.empty()) {
        std::cout << "\033[36mMitochondrion species contributing to osmotic pressure: ";
        for (size_t i = 0; i < includedMitoSpecies.size(); ++i) {
            std::cout << includedMitoSpecies[i];
            if (i < includedMitoSpecies.size() - 1) std::cout << ", ";
        }
        std::cout << "\033[0m\n";
    }

    return {osmoticPressureCytosol, osmoticPressureMito};
}

// Function to log all species and their concentrations with color coding
void logSpeciesWithColor(rr::RoadRunner& rr, std::map<std::string, double>& previousValues, bool logAllSpecies) {
    auto ids = rr.getFloatingSpeciesIds();
    auto concs = rr.getFloatingSpeciesConcentrationsV();

    // Access the ExecutableModel from RoadRunner
    auto* executableModel = rr.getModel();
    if (!executableModel) {
        throw std::runtime_error("No model loaded in RoadRunner.");
    }

    std::cout << "Floating species concentrations:\n";
    int columnCount = 0; // Track the number of columns printed
    for (size_t i = 0; i < ids.size(); ++i) {
        const std::string& id = ids[i];
        double currentValue = concs[i];
        double previousValue = previousValues[id];
        double diff = currentValue - previousValue;

        // Skip logging zero concentrations unless logAllSpecies is true
        if (!logAllSpecies && currentValue < 0.00004) {
            continue;
        }

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

        std::ostringstream entry;
        entry << id << " (" << compartment << "): "
            << std::fixed << std::setprecision(4) << currentValue
            << " (" << std::fixed << std::setprecision(4) << diff << ")";

        std::cout << std::left << std::setw(55) << entry.str() << "\033[0m";

        previousValues[id] = currentValue; // Update previous value

        // Print a newline after every 3 columns
        if (++columnCount % 3 == 0) {
            std::cout << "\n";
        }
    }
    if (columnCount % 3 != 0) {
        std::cout << "\n"; // Ensure the last line ends with a newline
    }

    // Log osmotic pressure with detailed calculation
    auto [osmoticPressureCytosol, osmoticPressureMito] = computeOsmoticPressure(rr);
    double R = 0.0821; // L·atm·mol⁻¹·K⁻¹
    double T = 298;    // Kelvin

    double c_cytosol = osmoticPressureCytosol / 1000; // mol/L
    double pressure_cytosol = c_cytosol * R * T;

    double c_mito = osmoticPressureMito / 1000; // mol/L
    double pressure_mito = c_mito * R * T;

    std::cout << "\033[34mOsmotic Pressure (Cytosol): " << std::fixed << std::setprecision(4)
            << c_cytosol << " mol/L * " << R << " * " << T << " = "
            << pressure_cytosol << " atm\033[0m\n";

    std::cout << "\033[36mOsmotic Pressure (Mitochondrion): " << std::fixed << std::setprecision(4)
            << c_mito << " mol/L * " << R << " * " << T << " = "
            << pressure_mito << " atm\033[0m\n";

    auto boundaryIds = rr.getBoundarySpeciesIds();
    std::cout << "Boundary species concentrations:\n";
    columnCount = 0; // Reset column count for boundary species
    for (const auto& id : boundaryIds) {
        double currentValue = rr.getValue(id);
        double previousValue = previousValues[id];
        double diff = currentValue - previousValue;

        // Skip logging zero concentrations unless logAllSpecies is true
        if (!logAllSpecies && currentValue == 0.0) {
            continue;
        }

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

        
        std::ostringstream entry;
        entry << id << " (" << compartment << "): "
            << std::fixed << std::setprecision(4) << currentValue;

        std::cout << std::left << std::setw(55) << entry.str() << "\033[0m";

        previousValues[id] = currentValue; // Update previous value

        // Print a newline after every 3 columns
        if (++columnCount % 3 == 0) {
            std::cout << "\n";
        }
    }
    if (columnCount % 3 != 0) {
        std::cout << "\n"; // Ensure the last line ends with a newline
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

    // Merge unit definitions from model2 into model1
    for (unsigned int i = 0; i < model2->getNumUnitDefinitions(); ++i) {
        auto unitDef = model2->getUnitDefinition(i)->clone();
        if (!model1->getUnitDefinition(unitDef->getIdAttribute())) {
            model1->addUnitDefinition(static_cast<libsbml::UnitDefinition*>(unitDef));
        } else {
            delete unitDef; // Avoid memory leak
        }
    }

    // Merge compartments from model2 into model1
    for (unsigned int i = 0; i < model2->getNumCompartments(); ++i) {
        auto compartment = model2->getCompartment(i)->clone();
        if (!model1->getCompartment(compartment->getIdAttribute())) {
            model1->addCompartment(static_cast<libsbml::Compartment*>(compartment));
        } else {
            delete compartment; // Avoid memory leak
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

    // Merge initial assignments from model2 into model1
    for (unsigned int i = 0; i < model2->getNumInitialAssignments(); ++i) {
        auto initAssign = model2->getInitialAssignment(i)->clone();
        if (!model1->getInitialAssignment(initAssign->getSymbol())) {
            model1->addInitialAssignment(static_cast<libsbml::InitialAssignment*>(initAssign));
        } else {
            delete initAssign; // Avoid memory leak
        }
    }

    // Merge rules from model2 into model1
    for (unsigned int i = 0; i < model2->getNumRules(); ++i) {
        auto rule = model2->getRule(i)->clone();
        if (!model1->getRule(rule->getVariable())) {
            model1->addRule(static_cast<libsbml::Rule*>(rule));
        } else {
            delete rule; // Avoid memory leak
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
    model->unsetMetaId();

    for (unsigned int i = 0; i < model->getNumCompartments(); ++i) {
        model->getCompartment(i)->unsetMetaId();
        model->getCompartment(i)->unsetAnnotation();
        model->getCompartment(i)->unsetNotes();
    }

    for (unsigned int i = 0; i < model->getNumSpecies(); ++i) {
        auto species = model->getSpecies(i);
        species->unsetMetaId();
        species->unsetAnnotation();
        species->unsetNotes();
        species->unsetMetaId();
        species->unsetSBOTerm();
    }

    for (unsigned int i = 0; i < model->getNumParameters(); ++i) {
        auto parameter = model->getParameter(i);
        parameter->unsetMetaId();
        parameter->unsetAnnotation();
        parameter->unsetNotes();
        parameter->unsetMetaId();
        parameter->unsetSBOTerm();
    }

    for (unsigned int i = 0; i < model->getNumReactions(); ++i) {
        auto reaction = model->getReaction(i);
        reaction->unsetMetaId();
        reaction->unsetAnnotation();
        reaction->unsetNotes();
        reaction->unsetMetaId();
        reaction->unsetSBOTerm();

        // Remove metaid from local parameters in the reaction's kinetic law
        if (auto kineticLaw = reaction->getKineticLaw()) {
            kineticLaw->unsetMetaId();
            for (unsigned int j = 0; j < kineticLaw->getNumParameters(); ++j) {
                kineticLaw->getParameter(j)->unsetMetaId();
            }
        }
        for(unsigned int j = 0; j < reaction->getNumReactants(); ++j) {
            auto reactant = reaction->getReactant(j);
            reactant->unsetMetaId();
            reactant->unsetAnnotation();
            reactant->unsetNotes();
        }
        for(unsigned int j = 0; j < reaction->getNumProducts(); ++j) {
            auto product = reaction->getProduct(j);
            product->unsetMetaId();
            product->unsetAnnotation();
            product->unsetNotes();
        }
        for(unsigned int j = 0; j < reaction->getNumModifiers(); ++j) {
            auto modifier = reaction->getModifier(j);
            modifier->unsetMetaId();
            modifier->unsetAnnotation();
            modifier->unsetNotes();
        }
    }

    for (unsigned int i = 0; i < model->getNumRules(); ++i) {
        model->getRule(i)->unsetMetaId();
        model->getRule(i)->unsetAnnotation();
        model->getRule(i)->unsetNotes();
    }

    for (unsigned int i = 0; i < model->getNumInitialAssignments(); ++i) {
        model->getInitialAssignment(i)->unsetMetaId();
        model->getInitialAssignment(i)->unsetAnnotation();
        model->getInitialAssignment(i)->unsetNotes();
    }

    for (unsigned int i = 0; i < model->getNumEvents(); ++i) {
        auto event = model->getEvent(i);
        event->unsetMetaId();
        event->unsetAnnotation();
        event->unsetNotes();
        event->unsetMetaId();
        event->unsetSBOTerm();
    }

    for (unsigned int i = 0; i < model->getNumUnitDefinitions(); ++i) {
        auto unitDef = model->getUnitDefinition(i);
        unitDef->unsetMetaId();
        unitDef->unsetAnnotation();
        unitDef->unsetNotes();
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

void checkDeathConditions(rr::RoadRunner& rr, bool& cellAlive, std::string& deathCause) {
    auto speciesIds = rr.getFloatingSpeciesIds();
    auto concentrations = rr.getFloatingSpeciesConcentrationsV();

    double NADH = -1.0;
    double NAD = -1.0;

    for (size_t i = 0; i < speciesIds.size(); ++i) {
        const std::string& id = speciesIds[i];
        double concentration = concentrations[i];

        if (id == "ATP" && concentration < 0.5) {
            deathCause = "Energy collapse";
            cellAlive = false;
            return;
        } 
        else if (id == "P" && concentration > 25.0) { // Pi overload
            deathCause = "Osmotic lysis";
            cellAlive = false;
            return;
        } 
        if (id == "ROS_m" && concentration > 0.2) { // threshold can be tuned
            deathCause = "Oxidative stress (ROS overload)";
            cellAlive = false;
            return;
        }
        else if (id == "NADH") {
            NADH = concentration;
        } 
        else if (id == "NAD") {
            NAD = concentration;
        }
    }

    // Only check Redox collapse after scanning all species
    if (NAD >= 0 && NADH >= 0) {
        double total = NAD + NADH;
        if (total > 0.0 && (NADH / total) > 0.9) { // Highly reduced state
            deathCause = "Redox collapse";
            cellAlive = false;
            return;
        }
    }
}

void logDebugInfo(rr::RoadRunner& rr, double t, double dt) {
    // Reaction fluxes
    auto reactionIds = rr.getReactionIds();
    auto reactionRates = rr.getReactionRates();
    std::cerr << "=== Reaction fluxes ===\n";
    for (size_t i = 0; i < reactionIds.size(); ++i) {
        std::cerr << std::setw(12) << reactionIds[i]
                  << ": " << std::fixed << std::setprecision(4)
                  << reactionRates[i] << "\n";
    }

    // Species rates-of-change
    auto speciesIds = rr.getFloatingSpeciesIds();
    auto speciesRates = rr.getRatesOfChange();
    auto speciesConcs = rr.getFloatingSpeciesConcentrationsV();
    std::cerr << "\n=== Species d[X]/dt ===\n";
    for (size_t i = 0; i < speciesIds.size(); ++i) {
        double newConcentration = speciesConcs[i] + speciesRates[i] * dt; // Estimate next concentration
        if (newConcentration < 0) {
            std::cerr << "\033[31m"; // Red color for negative concentration
        }
        std::cerr << std::setw(12) << speciesIds[i]
                  << ": d[X]/dt = " << std::fixed << std::setprecision(4) << speciesRates[i]
                  << ", [X] = " << speciesConcs[i] << "\033[0m\n"; // Reset color
    }

    // Full stoichiometry matrix
    ls::DoubleMatrix S = rr.getFullStoichiometryMatrix();
    std::cerr << "\n=== Biggest per-reaction contributor to each species ===\n";
    for (size_t i = 0; i < speciesIds.size(); ++i) {
        double maxContrib = 0.0;
        std::string culprit;
        for (size_t j = 0; j < reactionIds.size(); ++j) {
            double coeff = S(i, j);
            double contrib = coeff * reactionRates[j];
            if (std::abs(contrib) > std::abs(maxContrib)) {
                maxContrib = contrib;
                culprit = reactionIds[j];
            }
        }
        std::cerr << std::setw(12) << speciesIds[i]
                  << " <- " << std::setw(8) << culprit
                  << " (dXdt=" << maxContrib << ")\n";
    }

    // Largest flux and rate-of-change
    auto maxFluxIt = std::max_element(reactionRates.begin(), reactionRates.end(),
                                      [](double a, double b) { return std::abs(a) < std::abs(b); });
    size_t maxFluxIdx = std::distance(reactionRates.begin(), maxFluxIt);
    std::cerr << "\n>> Largest single flux: "
              << reactionIds[maxFluxIdx] << " = " << *maxFluxIt << "\n";

    auto maxRateIt = std::max_element(speciesRates.begin(), speciesRates.end(),
                                      [](double a, double b) { return std::abs(a) < std::abs(b); });
    size_t maxRateIdx = std::distance(speciesRates.begin(), maxRateIt);
    std::cerr << ">> Largest |dX/dt|: "
              << speciesIds[maxRateIdx] << " = " << *maxRateIt << "\n";
}

// Function to simulate lookahead and predict crash time
double simulateLookahead(rr::RoadRunner& rr, double currentTime, double lookaheadTime, double dt) {
    rr::RoadRunner rrLookahead = rr; // Clone the current RoadRunner instance
    double t = currentTime;
    bool cellAlive = true;
    std::string deathCause;

    while (t < currentTime + lookaheadTime) {
        try {
            rrLookahead.oneStep(t, dt);
            t += dt;

            // Check death conditions
            checkDeathConditions(rrLookahead, cellAlive, deathCause);
            if (!cellAlive) {
                std::cout << "->> Predicting death at t = " << t
                          << " (relative to current time: " << t - currentTime << ")"
                          << " due to: " << deathCause << std::endl;
                return t; // Return the predicted crash time
            }
        } catch (const std::exception& e) {
            std::cerr << "->> Predicting crash at t = " << t
                      << " (relative to current time: " << t - currentTime << ")"
                      << " with dt = " << dt << "\n"
                      << "Error: " << e.what() << "\n";
            return t;
        }
    }

    std::cout << "->> No crash predicted within the next " << lookaheadTime << " time units.\n";
    return -1.0; // Indicate no crash predicted
}

int main() {
    const std::string glyPath = "glycolysis.xml";
    const std::string tcaPath = "tca_cycle.xml";
    const std::string oxphosPath = "oxphos.xml";
    const std::string pppPath = "ppp.xml";

    if (!checkFileExists(glyPath) || !checkFileExists(tcaPath) || !checkFileExists(oxphosPath) || !checkFileExists(pppPath)) {
        return 1;
    }

    try {
        // Load and merge SBML models
        std::string sbml1 = loadSBMLFromFile(glyPath);
        std::string sbml2 = loadSBMLFromFile(tcaPath);
        std::string sbml3 = loadSBMLFromFile(oxphosPath);
        std::string sbml4 = loadSBMLFromFile(pppPath);
        std::string mergedSBML = mergeSBMLModels(mergeSBMLModels(mergeSBMLModels(sbml1, sbml2), sbml3), sbml4);

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

        auto integrator = rr.getIntegrator();
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

        integrator->setValue("relative_tolerance", 1e-6);
        integrator->setValue("absolute_tolerance", 1e-8);
        integrator->setValue("stiff", true);  // uses BDF
        // integrator->setValue("maximum_bdf_order", 5); // Optional but safe
        // integrator->setValue("maximum_adams_order", 12); // Optional
        // integrator->setValue("maximum_num_steps", 10000);  // default is 500
        integrator->setValue("minimum_time_step", 1e-10);

        std::cout << "Running real-time simulation. Press 'p' to play/pause, 'q' to quit.\n";

        double t = 0.0;
        double dt = 0.05;
        bool isPlaying = false; // Play/pause toggle
        bool logAllSpecies = true; // Flag to toggle between logging all species and non-zero species

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

        bool cellAlive = true;
        std::string deathCause;

        double lookaheadTime = 1.0; // Time units to look ahead
        double lookaheadDt = dt;  // Time step for lookahead simulation

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
                    std::cout << "Enter reaction ID to toggle: ";
                    std::string reactionId; // = "vAK";
                    std::cin >> reactionId;
                    toggleReactionRate(rr, reactionId);
                } else if (ch == 'v' || ch == 'V') {
                    logAllSpecies = !logAllSpecies; // Toggle logging mode
                    std::cout << "Logging mode: " << (logAllSpecies ? "All species" : "Non-zero species") << std::endl;
                } else if (ch == 'l') {
                    logDebugInfo(rr, t, dt);
                } else if (ch == 'f' || ch == 'F') {
                    // Trigger lookahead simulation
                    simulateLookahead(rr, t, lookaheadTime, lookaheadDt);
                }
            }

            if (isPlaying) {
                try {
                    
                    // Start performance timer
                    auto start = std::chrono::high_resolution_clock::now();
                    
                    rr.oneStep(t, dt);
                    
                    // Stop performance timer
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
                    
                    std::cout << "\033[2J\033[H"; // ANSI escape codes to clear screen and move cursor to top
                    // Log time taken for oneStep
                    std::cout << "Time taken for oneStep: " << (double)duration / 1000.0 << " milliseconds. Current time: " << t << "\n";
                    t += dt;
                } catch (const std::exception& e) {
                    std::cerr << "⚠️ CVODE failed at t = " << t
                              << " with dt = " << dt << "\n"
                              << "Error: " << e.what() << "\n\n";
                    logDebugInfo(rr, t, dt);
                    dt = std::max(dt * 0.5, 1e-6);
                    isPlaying = false;
                }

                // Clamp all floating species concentrations to zero if they go negative
                clampNegativeConcentrations(rr);

                // Clear the console and log species concentrations
                logSpeciesWithColor(rr, previousValues, logAllSpecies);

                // Run lookahead simulation every frame
                // simulateLookahead(rr, t, lookaheadTime, lookaheadDt);

                checkDeathConditions(rr, cellAlive, deathCause);

                if (!cellAlive) {
                    logDebugInfo(rr, t, dt);                    
                    std::cout << "Simulation ended. Death condition met: " << deathCause << std::endl;           
                    isPlaying = false;
                    break;
                }

                std::this_thread::sleep_for(std::chrono::milliseconds((long)(dt * 1000))); // Adjust for desired update rate
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

# Test project for libRoadRunner sim

## Overview

This project demonstrates real-time biochemical pathway simulation using the libRoadRunner library, a high-performance SBML (Systems Biology Markup Language) simulation engine. The application models cellular metabolism by integrating multiple interconnected metabolic pathways including glycolysis, the TCA cycle, oxidative phosphorylation, and the pentose phosphate pathway.

The simulator performs dynamic time-step integration of metabolic reactions, tracking the concentrations of dozens of metabolic species (glucose, ATP, NAD+, pyruvate, etc.) across cellular compartments (cytosol and mitochondrion). It calculates reaction fluxes, monitors energy states, and computes osmotic pressure differentials between compartments. The implementation showcases advanced features like SBML model merging, real-time keyboard-driven simulation control, and detection of metabolic failure conditions such as energy collapse.

This project serves as a practical example of using C++ with libRoadRunner for computational systems biology applications, demonstrating how to load SBML models, run simulations with custom integrators, query species concentrations, and analyze metabolic flux distributions in living systems.

## Setup

Grab the latest .zip **Release** from https://github.com/sys-bio/roadrunner/releases and extract it to `externals/libroadrunner-release`.

## Development

Start VSCode from the Developer Command Prompt for VS 2022

## Build and run

To just use CMakeLists, manually create a `build` directory and run the following commands from the root of the project:
```bash
mkdir build
cd build
cmake --build . --target run
```

Or to use CMakePresets.json, run the following command from root:

```bash
cmake --preset default
cmake --build --preset default
```

For the full command to build and run the project, use the following command:

```bash
cmake --build --preset default --config Release; .\build\Release\hello_roadrunner.exe
```

To run in Debug mode, remember to also change the CMakeLists
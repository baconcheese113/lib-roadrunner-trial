# Test project for libRoadRunner sim

## Setup

Grab the latest .zip **Release** from https://github.com/sys-bio/roadrunner/releases and extract it to `externals/libroadrunner-release`.

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
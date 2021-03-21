# Simple Variational Monte Carlo solve for FYS4411

### Compiling the project using CMake:
You can install CMake through one of the Linux package managers, e.g., `apt install cmake`, `pacman -S cmake`, etc. For Mac you can install using `brew install cmake`. Other ways of installing are shown here: [https://cmake.org/install/](https://cmake.org/install/).

In a Linux/Mac terminal this can be done by the following commands
```bash
# Create build-directory
mkdir build

# Move into the build-directory
cd build

# Run CMake to create a Makefile
cmake ../

# Make the Makefile using two threads
make -j2

# Move the executable to the top-directory
mv vmc ..
```
Or, simply run the script `compile_project` via
```bash
./compile_project
```
and the same set of commands are done for you. Now the project can be run by executing
```bash
./vmc
```
in the top-directory.

#### Cleaning the directory
Run `make clean` in the top-directory to remove the executable `vmc` and the `build`-directory.

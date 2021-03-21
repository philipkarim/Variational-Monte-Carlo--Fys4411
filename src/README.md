# Simple Variational Monte Carlo solve for FYS4411
## Code structure
### Main program
All parameters are defined in `main`. The most fundamental parts of the code are `system.cpp` which are executing the metropolis steps and `sampler.cpp` which is sampling the results. The codes in the `initialStates` folder are initializing the system by distributing the particles. The codes in the `wavefunction` folder are implementing exactly that, the wavefunction. The codes in the `Hamiltonian` folder is used to compute the hamiltonian energy.

### Analytics
The codes used for plotting the results are `plotting.py` which plots all the figures, except the one body denisty which is plotted by `onebodydensityplotter.py`. The last file used for analytics are `statisticalhandling.py` which uses the blocking method written by Marius Jonsson, and can be found [`here`](https://github.com/computative/block/blob/master/python/tictoc.py).

## Compiling the project using CMake:
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

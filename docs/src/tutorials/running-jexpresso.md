## Setup with CPUs

```bash
>> cd $JEXPRESSO_HOME
>> julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.API.precompile()"
```
followed by the following:

Push problem name to ARGS
You need to do this only when you run a new problem
```bash
push!(empty!(ARGS), EQUATIONS::String, EQUATIONS_CASE_NAME::String);
include("./src/Jexpresso.jl")
```

* PROBLEM_NAME is the name of your problem directory as $JEXPRESSO/problems/equations/problem_name
* PROBLEM_CASE_NAME is the name of the subdirectory containing the specific setup that you want to run: 

The path would look like 
```$JEXPRESSO/problems/equations/PROBLEM_NAME/PROBLEM_CASE_NAME```

Some examples available in this branch:

Example 1: to solve the 2D Euler equations with buyoancy and two passive tracers defined in [`problems/equations/CompEuler/thetaTracers`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/CompEuler/thetaTracers)  you would do the following:
```bash
push!(empty!(ARGS), "CompEuler", "thetaTracers");
include("./src/Jexpresso.jl")
```


Example 2: to solve the 3D Euler equations with buyoancy defined in [`problems/equations/CompEuler/3d`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/CompEuler/thetaTracers) you would do the following:
```bash
push!(empty!(ARGS), "CompEuler", "3d");
include("./src/Jexpresso.jl")
```


## Setup and Run with MPI

JEXPRESSO supports parallel execution using either OpenMPI or MPICH. Follow these steps to configure and run with your preferred MPI implementation.

### 1. Install MPI Implementation

Choose either OpenMPI or MPICH:

#### OpenMPI Installation
```bash
# Ubuntu/Debian
sudo apt install libopenmpi-dev openmpi-bin

# macOS (Homebrew)
brew install open-mpi

# Verify installation
mpiexec --version
```

#### MPICH Installation
```bash
# Ubuntu/Debian
sudo apt install mpich libmpich-dev

# macOS (Homebrew) 
brew install mpich

# Verify installation
mpiexec --version
```

### 2. Configure MPI Preferences

#### Automatic Configuration (Default Path)
Use this command when MPI (OpenMPI/MPICH) is installed in standard system paths (`/usr/bin`, `/usr/local/bin`, etc.):
```bash
julia --project=. -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; MPIPreferences.use_system_binary()'
```

#### Manual Configuration (For Multiple MPI Installations or MPI not in Default Path)
For MPI installations in non-standard locations (e.g., /opt/openmpi, or custom paths):
```bash
julia --project=. -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; MPIPreferences.use_system_binary(;extra_paths = ["/where/your/mpi/lib"])'
```
If MPI is installed via homebrew on macOS, the MPI lib path is:
```bash
/opt/homebrew/lib
```

### 3. Running with MPI

#### Basic Execution
```bash
mpiexec -n <NPROCS> julia --project=. -e 'push!(empty!(ARGS), "<EQUATIONS>", "<CASE_NAME>"); include("./src/Jexpresso.jl")'
```

#### Implementation-Specific Examples
```bash
mpiexec -n 4 julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "3d"); include("./src/Jexpresso.jl")'
```


### Troubleshooting

- **Library conflicts:** Clear existing preferences:
  ```bash
  rm -f LocalPreferences.toml
  ```
- **Path issues:** Verify paths with:
  ```bash
  which mpiexec
  which mpirun
  ```
  You may have to use the full aboslute path to mpiexec or mpirun and to julia like this if necessary:
  ```
  /opt/homebrew/Cellar/open-mpi/5.0.6/bin/mpirun -n 4 /Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "theta"); include("./src/Jexpresso.jl")'
  ```


- **Version mismatches:** Ensure consistent versions:
  ```bash
  mpicc --version
  mpif90 --version
  ```
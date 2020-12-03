# HPAI project

Date of project initialization: 2020-11-30

## Setup and installation
* Julia version 1.5.0
* Python 3.8.3
* This was run on a 16-inch MacBook Pro running MacOS Catalina 10.15.7

Activate the Conda env specified by the `environment.yml` file. Activate it in the Julia package `Conda.jl` by doing the following:
```julia
julia> ENV["CONDA_JL_HOME"] = "/path/to/miniconda/envs/hpai"  # change this to your path

pkg> build Conda
```

## How to run:
Using the enviroment above, execute the `main.jl` script file.

# ED_Sectors

ED_Sectors is a exact diagonalization code that is designed to be as painless as possible. Its main purpose is to diagonalize Hamiltonians for spin-1/2 chains, taking advantage of symmetries to speed up the computations. The code is optimized to have easy-to-use syntax and for convenience rather than speed; don't expect to get much beyond L = 20. 

## Installing 

ED_Sectors is written in Julia, a new language that combines python-like syntax with speed comperable to C and is designed for scientific computation. Eventually ED_Sectors might become a Julia package, but for now it's just a collection of files that can be included in any Julia program.

To install ED_Sectors:
1. Install Julia 0.7 or later from <https://julialang.org/downloads/>. Note that Julia is still in development and some bugs are still present. Linux repositories often have older versions of Julia, so it is advised that you use the latest version from the Julia website.
2. Clone this package by navigating to a directory of your choice and running

```
git clone https://github.com/DanielEricParker/ED_sectors
```

3. Include "main.jl" in your julia program via 

```
include("<directory of ED_sectors>/main.jl")
``` 
to get access to all ED_Sector commands.

## Getting Started

After you've installed ED_Sectors, check out the Examples to get started. They are well-documented and use most of the front-facing commands in the package. 

## Bugs
This is still a work-in-progress, and many bugs are probably lurking around still. Please start an issue if you find one! 


